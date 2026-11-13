function metrics = adaptive_link_metrics2(methodOverride, calOpts)
% Run ONE strategy per call. Set adaptationMethod below (or pass the method
% name as the first argument, e.g. adaptive_link_metrics2('oracle')).
%
% Adaptive methods : immediate, threshold, per_modulation, history, window,
%                    aimd, score, goodput_probe, q_learning
% Benchmarks       : always_bpsk   (most robust, non-adaptive)
%                    always_64qam  (highest rate, non-adaptive)
%                    oracle        (knows the instantaneous ACI state)
%                    always_qpsk, always_16qam (extra fixed modes)
%
% calOpts is used internally by the oracle for calibration; do not pass it.

    % --- System Parameters ---
    adaptationMethod = 'always_64qam'; % <-- CHANGE THIS to run a different method
    if nargin >= 1 && ~isempty(methodOverride)
        adaptationMethod = methodOverride;
    end
    calibrating = (nargin >= 2) && ~isempty(calOpts);

    rngSeed = 12345;   % same seed -> same ACI schedule for every method

    % Oracle: calibrate the best modulation for ACI OFF and ACI ON
    oracleMap = {};
    if strcmpi(adaptationMethod, 'oracle') && ~calibrating
        oracleMap = calibrate_oracle();
    end
    rng(rngSeed);

    modScheme = adaptive_controller([], adaptationMethod, true);    % 'BPSK', 'QPSK', '16-QAM', '64-QAM'
    numBlocks = 1000;
    numSymbols = 12;
    if calibrating
        numBlocks = calOpts.numBlocks;
    end

    % Baseband parameters
    Rs = 10e6;
    Bandwidth = 18e6;
    sps = 16;
    Fs = Rs * sps;

    %=========================================================
    %                    ACI PARAMETERS
    %=========================================================

    interference_type = 'ratioed';  % 'ratioed' or 'randomstr'

    % Mode 1: Ratioed periodic ACI schedule
    % [ACI_ON_blocks, ACI_OFF_blocks]
    aci_block_pattern = [10 20];

    % Mode 2: Random burst-based ACI schedule
    % [maximum duty cycle in percent, maximum consecutive ON blocks]
    randomstr_config = [40 30];

    max_power = 1;
    freq_separation = 10e6;

    % Amplifiers
    ampGain_Tx = 2;
    ampGain_Rx_linear = 1;
    ampGain_Rx_non_linear = 0.15;

    % Cubic nonlinearity
    a1 = 1;

    % AWGN variance
    noiseVar = 2.055e-21;

    %=========================================================
    %          VALIDATE ACI MODE AND GENERATE SCHEDULE
    %=========================================================

    aciSchedule = false(numBlocks, 1);

    if calibrating
        % Forced ACI state for oracle calibration
        aciSchedule(:) = calOpts.aciOn;
    else
        switch lower(interference_type)

            case 'ratioed'

                if numel(aci_block_pattern) ~= 2 || ...
                        any(aci_block_pattern < 0) || ...
                        any(aci_block_pattern ~= floor(aci_block_pattern)) || ...
                        sum(aci_block_pattern) == 0

                    error(['aci_block_pattern must be [nonnegative ACI_ON_blocks, ' ...
                           'nonnegative ACI_OFF_blocks], not both zero.']);
                end

                aci_on_blocks = aci_block_pattern(1);
                aci_off_blocks = aci_block_pattern(2);
                aci_cycle_len = aci_on_blocks + aci_off_blocks;

                aciSchedule = mod((1:numBlocks)' - 1, aci_cycle_len) < aci_on_blocks;

            case 'randomstr'

                if numel(randomstr_config) ~= 2 || ...
                        randomstr_config(1) < 0 || ...
                        randomstr_config(1) > 100 || ...
                        randomstr_config(1) ~= floor(randomstr_config(1)) || ...
                        randomstr_config(2) < 1 || ...
                        randomstr_config(2) ~= floor(randomstr_config(2))

                    error(['randomstr_config must be [duty_cycle_percent, ' ...
                           'positive_max_burst_length], with duty cycle from 0 to 100.']);
                end

                maxDutyCycle = randomstr_config(1);
                maxBurstLength = randomstr_config(2);

                maxInterferenceBlocks = floor(maxDutyCycle * numBlocks / 100);

                % Generate one random burst schedule for the complete run.
                bSch = 1;
                interferenceBlocksUsed = 0;

                while bSch <= numBlocks && ...
                        interferenceBlocksUsed < maxInterferenceBlocks

                    if rand < 0.5

                        remainingBlocks = numBlocks - bSch + 1;

                        remainingBudget = maxInterferenceBlocks - ...
                            interferenceBlocksUsed;

                        burstLength = randi([1, min([remainingBlocks, ...
                            remainingBudget, maxBurstLength])]);

                        aciSchedule(bSch:bSch + burstLength - 1) = true;

                        interferenceBlocksUsed = interferenceBlocksUsed ...
                            + burstLength;

                        % Insert a randomly sized clean gap after this burst.
                        cleanGap = randi([1, ...
                            max(1, numBlocks - bSch - burstLength + 1)]);

                        bSch = bSch + burstLength + cleanGap;

                    else

                        bSch = bSch + 1;

                    end
                end

            otherwise

                error('Unknown interference_type: %s. Use ratioed or randomstr.', ...
                    interference_type);

        end
    end

    % --- Modulation Setup ---
    switch modScheme
        case 'BPSK'
            M = 2;
        case 'QPSK'
            M = 4;
        case '16-QAM'
            M = 16;
        case '64-QAM'
            M = 64;
        otherwise
            error('Unsupported modulation scheme');
    end

    k = log2(M);
    blockSize = numSymbols * k;

    % --- RRC Filter Parameters ---
    span = 6;
    rolloff = (Bandwidth / Rs) - 1;

    hRRC = rcosdesign(rolloff, span, sps, 'sqrt');
    hRRC = hRRC / norm(hRRC);

    % (Live constellation and PSD figures are disabled to keep the run
    %  lightweight; the old commented-out plotting code was removed.)

    %=========================================================
    %                  SIMULATION LOOP
    %=========================================================

    nackHistory = false(numBlocks, 1);
    ackHistory = false(numBlocks, 1);
    aciHistory = false(numBlocks, 1);
    bitErrorsPerBlock = zeros(numBlocks, 1);
    bitsPerBlock = zeros(numBlocks, 1);
    blockTime_s = zeros(numBlocks, 1);
    usedModScheme = cell(numBlocks, 1);

    current_time = 0;

    for b = 1:numBlocks

        %% ---------------------------------------------------
        % 1. TRANSMITTER
        % ---------------------------------------------------

        % Oracle: choose modulation from the TRUE ACI state of this block.
        if strcmpi(adaptationMethod, 'oracle')
            modScheme = oracleMap{aciSchedule(b) + 1};   % {OFF, ON}
        end

        currentModScheme = modScheme;

        usedModScheme{b} = currentModScheme;
        completedBlocks = b;

        switch modScheme
            case 'BPSK'
                M = 2;
            case 'QPSK'
                M = 4;
            case '16-QAM'
                M = 16;
            case '64-QAM'
                M = 64;
        end

        k = log2(M);
        blockSize = numSymbols * k;

        txBits = randi([0 1], blockSize, 1);

        txSyms = encoder(txBits, modScheme, M);

        % RRC pulse shaping of the desired signal
        txWaveform = upfirdn(txSyms, hRRC, sps);
        txWaveform = txWaveform(:);

        txBB = txWaveform;

        desiredPower = mean(abs(txBB).^2);

        %----------------------------------------------------
        % Adjacent-channel interference generation
        %----------------------------------------------------

        interferencePower = max_power * desiredPower;

        % Select the ACI modulation independently for every block
        aciSchemes = {'BPSK', 'QPSK', '16-QAM', '64-QAM'};

        aciModIdx = randi(numel(aciSchemes));

        aciModScheme = aciSchemes{aciModIdx};

        switch aciModScheme
            case 'BPSK'
                Maci = 2;
            case 'QPSK'
                Maci = 4;
            case '16-QAM'
                Maci = 16;
            case '64-QAM'
                Maci = 64;
        end

        kAci = log2(Maci);

        aciBlockSize = numSymbols * kAci;

        aci_bits = randi([0 1], aciBlockSize, 1);

        aci_syms = encoder(aci_bits, aciModScheme, Maci);

        % RRC pulse shaping for ACI
        aci_waveform = upfirdn(aci_syms, hRRC, sps);
        aci_waveform = aci_waveform(:);

        L = min(length(txBB), length(aci_waveform));

        txBB = txBB(1:L);
        aci_waveform = aci_waveform(1:L);

        % Continuous time tracking
        t = current_time + (0:L-1)' / Fs;

        current_time = current_time + L / Fs;
        blockTime_s(b) = current_time;

        % Generate complex baseband ACI shifted in frequency
        aci_BB_raw = aci_waveform .* ...
            exp(1j * 2 * pi * freq_separation * t);

        rawAciPower = mean(abs(aci_BB_raw).^2);

        % Scale ACI to requested power
        aci_BB = sqrt(interferencePower / max(rawAciPower, eps)) ...
            .* aci_BB_raw;

        %----------------------------------------------------
        % ACI block schedule (ratioed / randomstr, built above)
        %----------------------------------------------------

        aciBlockOn = aciSchedule(b);

        aciHistory(b) = aciBlockOn;

        N = numel(aci_BB);

        if aciBlockOn
            aciMask = ones(N, 1);
            aciStatusText = 'ON';
        else
            aciMask = zeros(N, 1);
            aciStatusText = 'OFF';
        end

        aci_interference = aciMask .* aci_BB;

        % Actual transmitted baseband waveform includes ACI
        txBB = txBB + aci_interference;

        if ~calibrating
            fprintf('Block %d: ACI = %s | ACI modulation = %s | TX = %s\n', ...
                b, aciStatusText, aciModScheme, currentModScheme);
        end

        %% 2. TX AMPLIFIER
        txAmp = txBB * ampGain_Tx;

        %% 3. AWGN CHANNEL
        noise = sqrt(noiseVar/2) * ...
            (randn(size(txAmp)) + 1j*randn(size(txAmp)));

        channelGain = 0.5;

        rxWaveform = channelGain * txAmp + noise;

        %% 4. LINEAR RX AMPLIFIER
        rxAmp = rxWaveform * ampGain_Rx_linear;

        %% 5. NON-LINEAR RX AMPLIFIER
        % Baseband equivalent cubic nonlinearity: y = a1*x + a3*|x|^2*x
        rxNL = a1 * rxAmp + ...
            ampGain_Rx_non_linear * (abs(rxAmp).^2) .* rxAmp;

        rxWaveform = rxNL;

        %% 6. RECEIVER RRC FILTER
        rxFiltered = upfirdn(rxWaveform, hRRC, 1, sps);

        rxSyms = rxFiltered(span + 1 : span + numSymbols);

        %% 7. DEMODULATION
        rxBits = decoder(rxSyms, modScheme, M);

        %% 8. ACK / NACK
        [isAck, numErrors] = generate_ack_nack(txBits, rxBits);
        bitErrorsPerBlock(b) = numErrors;
        bitsPerBlock(b) = numel(txBits);

        ackHistory(b) = isAck;
        nackHistory(b) = ~isAck;

        % Modulation for the NEXT block (fixed for always_* methods,
        % overridden by the ACI state for 'oracle').
        modScheme = adaptive_controller(isAck, adaptationMethod, false);

    end

    %% =====================================================
    % FINAL PERFORMANCE METRICS
    % ======================================================

    schemeNames = {'BPSK', 'QPSK', '16-QAM', '64-QAM'};
    schemeRates = [1 2 4 6]; % Payload bits per modulation symbol

    % Bit error rate: total erroneous bits / total transmitted bits.
    totalBitErrors = sum(bitErrorsPerBlock(1:completedBlocks));
    totalTransmittedBits = sum(bitsPerBlock(1:completedBlocks));
    BER = totalBitErrors / max(totalTransmittedBits, 1);

    % Block error rate: fraction of blocks with one or more bit errors.
    blockErrorFlags = bitErrorsPerBlock(1:completedBlocks) > 0;
    BLER = mean(blockErrorFlags);
    ackRate = 1 - BLER;

    % Modulation selected for transmission in each block.
    schemeLevels = zeros(completedBlocks, 1);
    for b = 1:completedBlocks
        schemeLevels(b) = find(strcmp(schemeNames, ...
            usedModScheme{b}), 1);
    end

    % Fraction of transmitted blocks/time spent in each modulation mode.
    blockDurations = diff([0; blockTime_s(1:completedBlocks)]);
    totalSimulationTime = sum(blockDurations);
    modeFraction = zeros(1, numel(schemeNames));
    for mIdx = 1:numel(schemeNames)
        modeFraction(mIdx) = sum(blockDurations(schemeLevels == mIdx)) ...
            / totalSimulationTime;
    end

    % No retransmission of failed blocks in this simulator.
    actualRetransmissionRate = 0;
    nackRate = BLER; % Potential retransmission demand, if ARQ were enabled.

    % Useful payload delivered by ACKed blocks.
    blockRates = schemeRates(schemeLevels).';
    successfulBitsPerBlock = blockRates .* numSymbols .* ...
        double(~blockErrorFlags);
    totalSuccessfulBits = sum(successfulBitsPerBlock);

    usefulThroughput_bps = totalSuccessfulBits / totalSimulationTime;
    goodput_bitsPerSymbol = totalSuccessfulBits / ...
        (completedBlocks * numSymbols);

    time_s = blockTime_s(1:completedBlocks);

    % Return all scalar and per-block metrics.
    metrics = struct();
    metrics.Method = adaptationMethod;
    metrics.BER = BER;
    metrics.BLER = BLER;
    metrics.ActualRetransmissionRate = actualRetransmissionRate;
    metrics.NACKRate_RetransmissionDemand = nackRate;
    metrics.Goodput_bitsPerSymbol = goodput_bitsPerSymbol;
    metrics.UsefulThroughput_bps = usefulThroughput_bps;
    metrics.TotalBitErrors = totalBitErrors;
    metrics.TotalTransmittedBits = totalTransmittedBits;
    metrics.TotalSuccessfulBits = totalSuccessfulBits;
    metrics.ModeNames = schemeNames;
    metrics.ModeFractions = modeFraction;
    metrics.Time_s = time_s;
    metrics.SelectedModulation = usedModScheme(1:completedBlocks);
    metrics.BitErrorsPerBlock = bitErrorsPerBlock(1:completedBlocks);
    metrics.BlockErrorFlags = blockErrorFlags;
    metrics.ACI_Active = aciHistory(1:completedBlocks);
    if ~isempty(oracleMap)
        metrics.OraclePolicy = struct('ACI_OFF', oracleMap{1}, ...
                                      'ACI_ON',  oracleMap{2});
    end

    % Calibration runs only need the metrics: skip files, plots and prints.
    if calibrating
        return;
    end

    % Files are tagged with the method name so separate runs do not
    % overwrite each other.
    tag = adaptationMethod;

    summaryTable = table(BER, BLER, actualRetransmissionRate, nackRate, ...
        goodput_bitsPerSymbol, usefulThroughput_bps, ...
        modeFraction(1), modeFraction(2), modeFraction(3), modeFraction(4), ...
        'VariableNames', {'BER','BLER','Actual_Retransmission_Rate', ...
        'NACK_Rate_Retransmission_Demand','Goodput_bits_per_symbol', ...
        'Useful_Throughput_bps','Fraction_BPSK','Fraction_QPSK', ...
        'Fraction_16QAM','Fraction_64QAM'});
    writetable(summaryTable, ['adaptive_link_summary_' tag '.csv']);

    perBlockTable = table((1:completedBlocks).', time_s, ...
        string(usedModScheme(1:completedBlocks)), ...
        bitsPerBlock(1:completedBlocks), ...
        bitErrorsPerBlock(1:completedBlocks), ...
        blockErrorFlags, double(aciHistory(1:completedBlocks)), ...
        'VariableNames', {'Block','Time_s','Selected_Modulation', ...
        'Transmitted_Bits','Bit_Errors','Block_Error','ACI_Active'});
    writetable(perBlockTable, ['adaptive_link_per_block_' tag '.csv']);
    save(['adaptive_link_metrics_' tag '.mat'], 'metrics', 'summaryTable', ...
        'perBlockTable', 'adaptationMethod');

    % Compact dashboard with the requested outputs.
    fig = figure('Name',['Adaptive Link Performance Metrics - ' tag], ...
        'NumberTitle','off','Color','w','Visible','on');
    tiledlayout(fig, 2, 2, 'TileSpacing','compact','Padding','compact');

    % 1. Selected modulation versus time.
    nexttile;
    stairs(time_s * 1e3, schemeLevels, 'LineWidth', 1.4);
    grid on;
    yticks(1:4);
    yticklabels(schemeNames);
    xlabel('Time (ms)');
    ylabel('Selected modulation');
    title('Modulation versus Time');

    % 2. Fraction of time in each modulation.
    nexttile;
    bar(100 * modeFraction);
    grid on;
    xticks(1:4);
    xticklabels(schemeNames);
    ylabel('Time fraction (%)');
    ylim([0 100]);
    title('Time in Each Modulation Mode');

    % 3. BER, BLER and retransmission-related rates.
    nexttile;
    bar([BER, BLER, actualRetransmissionRate, nackRate]);
    grid on;
    xticks(1:4);
    xticklabels({'BER','BLER','Actual retx.','NACK demand'});
    ylabel('Rate');
    title('Error and Retransmission Metrics');

    % 4. Useful throughput and normalized goodput.
    nexttile;
    yyaxis left;
    bar(1, usefulThroughput_bps / 1e6);
    ylabel('Useful throughput (Mbit/s)');
    yyaxis right;
    plot(1, goodput_bitsPerSymbol, 'ro', ...
        'MarkerFaceColor','r','MarkerSize',7);
    ylabel('Goodput (bits/symbol)');
    xlim([0.5 1.5]);
    xticks(1);
    xticklabels({'Delivered payload'});
    grid on;
    title('Useful Throughput / Goodput');

    exportgraphics(fig, ['adaptive_link_metrics_' tag '.png'], 'Resolution', 160);

    fprintf('\n========== ADAPTIVE LINK METRICS ==========\n');
    fprintf('Adaptation method                 : %s\n', adaptationMethod);
    if ~isempty(oracleMap)
        fprintf('Oracle policy                     : ACI OFF -> %s | ACI ON -> %s\n', ...
            oracleMap{1}, oracleMap{2});
    end
    fprintf('Completed blocks                  : %d\n', completedBlocks);
    fprintf('ACK rate                          : %.2f%%\n', 100*ackRate);
    fprintf('Total successful payload          : %d bits\n', totalSuccessfulBits);
    fprintf('BER                               : %.6g\n', BER);
    fprintf('BLER                              : %.4f (%.2f%%)\n', BLER, 100*BLER);
    fprintf('Actual retransmission rate        : %.4f (no ARQ in current code)\n', ...
        actualRetransmissionRate);
    fprintf('NACK rate / retransmission demand : %.4f (%.2f%%)\n', ...
        nackRate, 100*nackRate);
    fprintf('Goodput                           : %.4f bits/symbol\n', ...
        goodput_bitsPerSymbol);
    fprintf('Useful throughput                 : %.4f Mbit/s\n', ...
        usefulThroughput_bps/1e6);
    fprintf('Time in BPSK / QPSK / 16QAM / 64QAM: %.1f%% / %.1f%% / %.1f%% / %.1f%%\n', ...
        100*modeFraction(1), 100*modeFraction(2), ...
        100*modeFraction(3), 100*modeFraction(4));
    fprintf('Saved CSV, MAT and PNG metric outputs in: %s\n', pwd);
    fprintf('===========================================\n');

end


%% =========================================================
%        ORACLE CALIBRATION (best modulation per ACI state)
% ==========================================================

function oracleMap = calibrate_oracle()
% Runs every fixed modulation with ACI forced OFF and then forced ON and
% returns {best modulation when OFF, best modulation when ON} by goodput.

    fixedMethods = {'always_bpsk','always_qpsk','always_16qam','always_64qam'};
    schemeNames  = {'BPSK','QPSK','16-QAM','64-QAM'};
    gp = zeros(2, 4);   % row 1 = ACI OFF, row 2 = ACI ON

    opts.numBlocks = 300;
    for s = 1:2
        opts.aciOn = (s == 2);
        for m = 1:4
            mt = adaptive_link_metrics2(fixedMethods{m}, opts);
            gp(s, m) = mt.Goodput_bitsPerSymbol;
        end
    end

    [~, bestOff] = max(gp(1, :));
    [~, bestOn]  = max(gp(2, :));
    oracleMap = {schemeNames{bestOff}, schemeNames{bestOn}};
end


%% =========================================================
%                    ENCODER
% ==========================================================

function syms = encoder(bits, modScheme, M)

    k = log2(M);

    bitMatrix = reshape(bits, k, [])';

    intData = bi2de(bitMatrix, 'left-msb');

    switch modScheme

        case 'BPSK'
            syms = pskmod(intData, 2);

        case 'QPSK'
            syms = pskmod(intData, 4, pi/4);

        otherwise
            syms = qammod(intData, M, ...
                'UnitAveragePower', true);

    end

end


%% =========================================================
%                    DECODER
% ==========================================================

function bits = decoder(syms, modScheme, M)

    switch modScheme

        case 'BPSK'
            rxInt = pskdemod(syms, 2);

        case 'QPSK'
            rxInt = pskdemod(syms, 4, pi/4);

        otherwise
            rxInt = qamdemod(syms, M, ...
                'UnitAveragePower', true);

    end

    k = log2(M);

    bitMatrix = de2bi(rxInt, k, 'left-msb');

    bits = reshape(bitMatrix', [], 1);

end


%% =========================================================
%                    ACK / NACK
% ==========================================================

function [isAck, errors] = generate_ack_nack(txBits, rxBits)

    errors = sum(txBits ~= rxBits);

    if errors == 0
        isAck = true;
    else
        isAck = false;
    end

end


%% =========================================================
%                ADAPTIVE MODULATION CONTROLLER
% ==========================================================

function modScheme = adaptive_controller(ack, method, reset)
% Adaptive modulation using ACK/NACK feedback.
% Methods: immediate, threshold, per_modulation, history, window, aimd,
%          score, goodput_probe, q_learning
% Benchmarks: always_bpsk, always_qpsk, always_16qam, always_64qam, oracle
%             (for 'oracle' the main code overrides the choice using the
%              true ACI state; the controller itself returns 64-QAM)
%
% ack    = 1 for ACK, 0 for NACK
% reset  = true to reset all controller memory
%
% Returns the modulation for the NEXT block.

persistent idx ackCount nackCount
persistent histAck
persistent blockCount lastSwitch
persistent aimdScore

persistent modHistory modObs probeCycle historyLocked
persistent bestModIdx
persistent holdRemaining
persistent scoreValue

% Goodput-aware probing state
persistent gpHistory gpBlocksSinceProbe
persistent gpProbeActive gpProbeBase gpProbeCount gpProbeAcks

% Algorithm 9: Q-Learning state
persistent qTable qCount qSince qPrevState
schemes = {'BPSK','QPSK','16-QAM','64-QAM'};

% Relative useful payload rates (bits/symbol)
rates = [1 2 4 6];

if nargin < 2 || isempty(method)
    method = 'threshold';
end
if nargin < 3
    reset = false;
end

% Algorithm 9: Online tabular Q-learning (no model, no pre-training)
qGamma    = 0.60;   % Discount (actions affect the next state via used mod)
qAlphaMin = 0.10;   % Learning-rate floor
qEps0     = 0.15;   % Initial exploration probability
qEpsMin   = 0.02;   % Minimum exploration probability
qEpsDecay = 0.99;   % Exploration decay per block
qMaxSince = 3;      % Cap on "blocks since last NACK"

% Initialize/reset
if isempty(idx) || reset
    idx = 4;
    if strcmpi(method, 'goodput_probe')
        idx = 1;  % Start from BPSK for probing / Q-learning
    end
    % Fixed (non-adaptive) benchmarks
    switch lower(method)
        case 'always_bpsk',  idx = 1;
        case 'always_qpsk',  idx = 2;
        case 'always_16qam', idx = 3;
        case 'always_64qam', idx = 4;
    end
    ackCount = 0;
    nackCount = 0;

    histAck = [];

    blockCount = 0;
    lastSwitch = 0;
    aimdScore = 1;
    modHistory = cell(1,4);
    modObs = zeros(1,4);

    for k = 1:4
        modHistory{k} = zeros(1,0);  % Initially no observations
    end
    bestModIdx = 3;

    probeCycle = 0;       % 1 -> 16-QAM, 2 -> QPSK, 3 -> BPSK
    historyLocked = false;
    holdRemaining = 0;

    % Algorithm 7 persistent state
    scoreValue = 5;

    % Algorithm 8 persistent state
    gpHistory = cell(1,4);
    for k = 1:4
        gpHistory{k} = [];
    end
    gpBlocksSinceProbe = 0;
    gpProbeActive = false;
    gpProbeBase = 1;
    gpProbeCount = 0;
    gpProbeAcks = 0;

    % Algorithm 9 state
    qTable = zeros(4*(qMaxSince+1), 4);
    for a = 1:4
        qTable(:,a) = rates(a)/(1-qGamma);   % Optimistic initialization
    end
    qCount = zeros(size(qTable));
    qSince = qMaxSince;                       % No recent NACK
    qPrevState = qMaxSince*4 + 4;             % (since = max, last mod = 64-QAM)
    if strcmpi(method, 'q_learning')
        idx = 4;
    end
end

% General parameters
m = 3;
n = 1;

% Algorithm 3: thresholds for BPSK, QPSK, 16-QAM, 64-QAM
mLevels = [3 3 4 5];
nLevels = [1 1 1 1];

% Algorithm 4: history-based inference and learning
historyWindow = 5;   % Observations required per modulation
holdDuration = 5;    % Blocks to use selected modulation after 64-QAM NACK


% Algorithm 5: sliding-window hysteresis
Wwindow = 5;
TL = 0.40;
TH = 0.80;
dwellTime = 2;

% Algorithm 6: AIMD
alpha = 0.5;
beta = 0.70;
xmin = 1;
xmax = 6;
levels = [1 2 4 6];

% Algorithm 7: Recent-failure score controller
scoreA = 0.5;
scoreB = 1.0;
scoreMin = 0;
scoreMax = 10;
scoreQL = 3;
scoreQH = 7;

% Algorithm 8: Goodput-aware probing
probeInterval = 10;
probeDuration = 5;
probeEpsilon = 0.20;
probeDelta = 0.10;
gpHistoryWindow = 20;

if isempty(ack)
    if nargin >= 3 && reset
        modScheme = schemes{idx};
        return;
    else
        error('ack must be either 0 (NACK) or 1 (ACK).');
    end
elseif ~isscalar(ack) || ~ismember(ack, [0 1])
    error('ack must be either 0 (NACK) or 1 (ACK).');
end

% Save the modulation used for the current block.
usedIdx = idx;
blockCount = blockCount + 1;


switch lower(method)

    %--------------------------------------------------------------
    % Benchmarks: non-adaptive fixed modulation (idx never changes)
    %--------------------------------------------------------------
    case {'always_bpsk','always_qpsk','always_16qam','always_64qam'}
        % Nothing to adapt: idx was fixed at reset.

    %--------------------------------------------------------------
    % Benchmark: oracle (choice is made in the main code from the
    % true ACI state; nothing to do here).
    %--------------------------------------------------------------
    case 'oracle'
        idx = 4;

    %--------------------------------------------------------------
    % Algorithm 1: Immediate ACK/NACK adaptation
    %--------------------------------------------------------------
    case 'immediate'
        if ack == 1
            idx = min(idx + 1, 4);
        else
            idx = max(idx - 1, 1);
        end

    %--------------------------------------------------------------
    % Algorithm 2: Consecutive m-ACK / n-NACK thresholds
    %--------------------------------------------------------------
    case 'threshold'
        if ack == 1
            ackCount = ackCount + 1;
            nackCount = 0;
        else
            nackCount = nackCount + 1;
            ackCount = 0;
        end

        if ackCount >= m
            idx = min(idx + 1, 4);
            ackCount = 0;
            nackCount = 0;
        elseif nackCount >= n
            idx = max(idx - 1, 1);
            ackCount = 0;
            nackCount = 0;
        end

    %--------------------------------------------------------------
    % Algorithm 3: Per-modulation thresholds
    %--------------------------------------------------------------
    case 'per_modulation'
        if ack == 1
            ackCount = ackCount + 1;
            nackCount = 0;
        else
            nackCount = nackCount + 1;
            ackCount = 0;
        end

        if ackCount >= mLevels(usedIdx)
            idx = min(idx + 1, 4);
            ackCount = 0;
            nackCount = 0;
        elseif nackCount >= nLevels(usedIdx)
            idx = max(idx - 1, 1);
            ackCount = 0;
            nackCount = 0;
        end


%--------------------------------------------------------------
% Algorithm 4: Cyclic exploration, then fixed best-mod selection
%--------------------------------------------------------------
case 'history'

    % Record ACK/NACK for the modulation used in this block.
    modHistory{usedIdx} = [modHistory{usedIdx}, ack];

    % Retain only the most recent historyWindow observations.
    if numel(modHistory{usedIdx}) > historyWindow
        modHistory{usedIdx} = ...
            modHistory{usedIdx}(end-historyWindow+1:end);
    end

    modObs(usedIdx) = min(modObs(usedIdx) + 1, historyWindow);

    %----------------------------------------------------------
    % Lock once every modulation has enough observations.
    % Select the best LOWER modulation only once.
    %----------------------------------------------------------
    if ~historyLocked && all(modObs >= historyWindow)

        candidates = [3 2 1];  % 16-QAM, QPSK, BPSK
        goodput = zeros(1, numel(candidates));

        for k = 1:numel(candidates)
            s = candidates(k);
            ackRatio = mean(modHistory{s});
            goodput(k) = rates(s) * ackRatio;
        end

        [~, bestPos] = max(goodput);
        bestModIdx = candidates(bestPos);

        historyLocked = true;

        fprintf(['History locked at block %d. Best modulation: %s, ' ...
                 'estimated goodput: %.3f\n'], ...
                 blockCount, schemes{bestModIdx}, goodput(bestPos));
    end

    %----------------------------------------------------------
    % Phase 1: Cyclic exploration (UNCHANGED)
    %----------------------------------------------------------
    if ~historyLocked

        if usedIdx == 4

            if ack == 0
                % A 64-QAM NACK triggers the next cyclic probe.
                probeCycle = mod(probeCycle, 3) + 1;

                % Cycle: 16-QAM -> QPSK -> BPSK -> repeat.
                probeSchemes = [3 2 1];
                idx = probeSchemes(probeCycle);
            else
                % 64-QAM ACK: continue with 64-QAM.
                idx = 4;
            end

        else
            % After a probe, always return to 64-QAM.
            idx = 4;
        end

    %----------------------------------------------------------
    % Phase 2: Locked operation using fixed best modulation
    %----------------------------------------------------------
    else

        if usedIdx == 4 && ack == 0

            % 64-QAM NACK: use the fixed selected modulation.
            % Do not recalculate goodput or inspect history.
            idx = bestModIdx;

            % The immediate next block counts as block 1.
            holdRemaining = holdDuration - 1;

        elseif usedIdx ~= 4 && holdRemaining > 0

            % Continue using the same selected modulation.
            holdRemaining = holdRemaining - 1;

            if holdRemaining == 0
                idx = 4;
            else
                idx = bestModIdx;
            end

        else

            % Normal operation: use 64-QAM.
            idx = 4;
            holdRemaining = 0;

        end
    end

%--------------------------------------------------------------
    % Algorithm 5: Sliding-window ACK ratio with hysteresis
    %--------------------------------------------------------------
    case 'window'

        % Record feedback for the current window.
        histAck = [histAck ack];

        % Keep only the most recent Wwindow observations.
        if numel(histAck) > Wwindow
            histAck = histAck(end-Wwindow+1:end);
        end

        % Wait until the window is completely filled.
        if numel(histAck) < Wwindow
            % Not enough samples: retain the current modulation.
            % Do not make an adaptation decision.

        else
            % Calculate ACK ratio over the full window.
            pACK = mean(histAck);

            % Check dwell time.
            canSwitch = (blockCount-lastSwitch) >= dwellTime;

            newIdx = idx;

            if canSwitch && pACK >= TH
                newIdx = min(idx+1, 4);

            elseif canSwitch && pACK <= TL
                newIdx = max(idx-1, 1);
            end

            % Reset the window only after an actual modulation change.
            if newIdx ~= idx
                idx = newIdx;
                lastSwitch = blockCount;
                histAck = [];
            end
        end
    %--------------------------------------------------------------
    % Algorithm 6: Additive-Increase / Multiplicative-Decrease
    %--------------------------------------------------------------
    case 'aimd'
        if ack == 1
            aimdScore = min(aimdScore+alpha,xmax);
        else
            aimdScore = max(beta*aimdScore,xmin);
        end

        % Map the score to the nearest supported modulation level.
        [~,idx] = min(abs(levels-aimdScore));

    %--------------------------------------------------------------
    % Algorithm 7: Recent-Failure Penalty / Score Controller
    %--------------------------------------------------------------
    case 'score'

        if ack == 1
            scoreValue = min(scoreValue + scoreA, scoreMax);
        else
            scoreValue = max(scoreValue - scoreB, scoreMin);
        end

        if scoreValue >= scoreQH
            idx = min(usedIdx + 1, 4);
        elseif scoreValue <= scoreQL
            idx = max(usedIdx - 1, 1);
        else
            idx = usedIdx;
        end
    %--------------------------------------------------------------
    % Algorithm 8: Goodput-aware probing
    %--------------------------------------------------------------
    case 'goodput_probe'

        % Store feedback for the modulation used in this block.
        gpHistory{usedIdx} = [gpHistory{usedIdx}, ack];

        if numel(gpHistory{usedIdx}) > gpHistoryWindow
            gpHistory{usedIdx} = ...
                gpHistory{usedIdx}(end-gpHistoryWindow+1:end);
        end

        % If a probe is active, collect its ACK/NACK observations.
        if gpProbeActive && usedIdx == gpProbeBase + 1

            gpProbeCount = gpProbeCount + 1;
            gpProbeAcks = gpProbeAcks + ack;

            if gpProbeCount >= probeDuration

                % Estimated probe performance.
                probeAckRatio = gpProbeAcks / gpProbeCount;
                probeFailureRate = 1 - probeAckRatio;
                probeGoodput = rates(usedIdx) * probeAckRatio;

                % Estimated baseline performance.
                baseHistory = gpHistory{gpProbeBase};

                if isempty(baseHistory)
                    baseAckRatio = 0;
                else
                    baseAckRatio = mean(baseHistory);
                end

                baseGoodput = rates(gpProbeBase) * baseAckRatio;

                % Accept only if goodput improves enough and reliability
                % satisfies the failure-rate target.
                if probeGoodput >= baseGoodput + probeDelta && ...
                        probeFailureRate <= probeEpsilon

                    idx = usedIdx;
                    fprintf(['Probe accepted: %s, goodput %.3f ' ...
                        '(baseline %.3f)\n'], ...
                        schemes{usedIdx}, probeGoodput, baseGoodput);

                else
                    idx = gpProbeBase;
                    fprintf(['Probe rejected: %s, goodput %.3f ' ...
                        '(baseline %.3f), failure rate %.1f%%\n'], ...
                        schemes{usedIdx}, probeGoodput, baseGoodput, ...
                        100 * probeFailureRate);
                end

                gpProbeActive = false;
                gpBlocksSinceProbe = 0;
                gpProbeCount = 0;
                gpProbeAcks = 0;
            else
                % Continue probing the same higher modulation.
                idx = usedIdx;
            end

        elseif gpProbeActive

            % Safety fallback if the current modulation differs from the
            % expected probe modulation.
            idx = gpProbeBase;
            gpProbeActive = false;
            gpBlocksSinceProbe = 0;

        else
            gpBlocksSinceProbe = gpBlocksSinceProbe + 1;
            idx = usedIdx;

            % Probe only when there is a higher modulation available and
            % a baseline estimate has been collected.
            if gpBlocksSinceProbe >= probeInterval && usedIdx < 4 && ...
                    ~isempty(gpHistory{usedIdx})

                gpProbeBase = usedIdx;
                gpProbeActive = true;
                gpProbeCount = 0;
                gpProbeAcks = 0;

                % The next block starts the probe.
                idx = usedIdx + 1;
            end
        end
    %--------------------------------------------------------------
    % Algorithm 9: Model warm-started belief-state Q-Learning
    %--------------------------------------------------------------
    case 'q_learning'
        % Update the "blocks since last NACK" counter
        if ack == 0
            qSince = 0;
        else
            qSince = min(qSince + 1, qMaxSince);
        end

        % Next state = (blocks since NACK, modulation just used)
        nextState = qSince*4 + usedIdx;

        % Reward = delivered payload bits/symbol
        reward = rates(usedIdx) * ack;

        % Q-update for the transition (qPrevState, usedIdx) -> nextState
        qCount(qPrevState, usedIdx) = qCount(qPrevState, usedIdx) + 1;
        alphaQ = max(qAlphaMin, 1/(1 + qCount(qPrevState, usedIdx)));
        target = reward + qGamma * max(qTable(nextState,:));
        qTable(qPrevState, usedIdx) = qTable(qPrevState, usedIdx) + ...
            alphaQ * (target - qTable(qPrevState, usedIdx));

        % Decaying epsilon-greedy action selection
        qEps = max(qEpsMin, qEps0 * qEpsDecay^blockCount);
        if rand < qEps
            idx = randi(4);
        else
            qValues = qTable(nextState,:);
            best = find(qValues >= max(qValues) - 1e-9);
            idx = best(randi(numel(best)));
        end

        qPrevState = nextState;
    otherwise
        error(['Unknown method. Use immediate, threshold, per_modulation, ' ...
               'history, window, aimd, score, goodput_probe, q_learning, ' ...
               'always_bpsk, always_qpsk, always_16qam, always_64qam, or oracle.']);
end


modScheme = schemes{idx};
end