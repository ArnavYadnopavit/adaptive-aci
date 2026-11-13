function sweep_snr_visible_plots()
% Sweep AWGN SNR and compare adaptive modulation algorithms.
% SNR is defined at RX input relative to desired signal power, before ACI.

    clc;

    % Force MATLAB figures to be visible (overrides hidden default settings).
    set(groot, 'defaultFigureVisible', 'on');

    snrValues = -5:1:30;             % dB
    methods = {'history', 'q_learning', 'immediate'};
    fSep = 10e6;                     % Fixed ACI frequency separation
    randomSeed = 12345;

    nSNR = numel(snrValues);
    nMethods = numel(methods);
    goodput = zeros(nSNR, nMethods);
    ackRate = zeros(nSNR, nMethods);
    successfulBits = zeros(nSNR, nMethods);

    fprintf('===== SNR SWEEP =====\n');
    for i = 1:nSNR
        for j = 1:nMethods
            rng(randomSeed, 'twister');
            evalc('[g, a, sb] = simulate_adaptive_link_snr(fSep, methods{j}, snrValues(i));');
            goodput(i,j) = g;
            ackRate(i,j) = a;
            successfulBits(i,j) = sb;
            fprintf('SNR = %5.1f dB | %-10s | Goodput = %.4f | ACK = %.2f%%\n', ...
                snrValues(i), methods{j}, g, 100*a);
        end
    end

    results = table(snrValues(:), ...
        goodput(:,1), goodput(:,2), goodput(:,3), ...
        100*ackRate(:,1), 100*ackRate(:,2), 100*ackRate(:,3), ...
        'VariableNames', {'SNR_dB', ...
        'Goodput_History', 'Goodput_QLearning', 'Goodput_Immediate', ...
        'ACK_History_pct', 'ACK_QLearning_pct', 'ACK_Immediate_pct'});

    writetable(results, 'snr_comparison.csv');
    save('snr_comparison.mat', 'results', 'snrValues', 'methods', ...
        'goodput', 'ackRate', 'successfulBits', 'fSep');

    % Only two lightweight plots, generated after all simulations.
    % Save PNG copies as a fallback if the MATLAB figure windows are hidden.
    fig1 = figure('Name', 'Goodput vs SNR', ...
        'NumberTitle', 'off', 'Color', 'w', 'Visible', 'on');
    plot(snrValues, goodput(:,1), '-o', 'LineWidth', 1.4); hold on;
    plot(snrValues, goodput(:,2), '-s', 'LineWidth', 1.4);
    plot(snrValues, goodput(:,3), '-^', 'LineWidth', 1.4);
    grid on;
    xlabel('SNR (dB)');
    ylabel('Average goodput (bits/symbol)');
    title('Adaptive Modulation Goodput vs SNR');
    legend(methods, 'Location', 'best');
    drawnow;
    exportgraphics(fig1, fullfile(pwd, 'goodput_vs_snr.png'), 'Resolution', 150);

    fig2 = figure('Name', 'ACK Rate vs SNR', ...
        'NumberTitle', 'off', 'Color', 'w', 'Visible', 'on');
    plot(snrValues, 100*ackRate(:,1), '-o', 'LineWidth', 1.4); hold on;
    plot(snrValues, 100*ackRate(:,2), '-s', 'LineWidth', 1.4);
    plot(snrValues, 100*ackRate(:,3), '-^', 'LineWidth', 1.4);
    grid on;
    xlabel('SNR (dB)');
    ylabel('ACK rate (%)');
    ylim([0 100]);
    title('ACK Rate vs SNR');
    legend(methods, 'Location', 'best');
    drawnow;
    exportgraphics(fig2, fullfile(pwd, 'ack_rate_vs_snr.png'), 'Resolution', 150);

    fprintf('\nSaved files in: %s\n', pwd);
    fprintf('  snr_comparison.csv\n');
    fprintf('  snr_comparison.mat\n');
    fprintf('  goodput_vs_snr.png\n');
    fprintf('  ack_rate_vs_snr.png\n');
end


function [goodput, ackRate, successfulBits] = ...
    simulate_adaptive_link_snr(fSep, method, snrDb)

    if nargin < 1 || isempty(fSep)
        fSep = 10e6;
    end

    if nargin < 2 || isempty(method)
        method = 'q_learning';
    end

    if nargin < 3 || isempty(snrDb)
        snrDb = 10;
    end

    if ~isscalar(snrDb) || ~isfinite(snrDb)
        error('SNR must be a finite scalar in dB.');
    end

    freq_separation = fSep;

    if ~isscalar(freq_separation) || ...
            ~isfinite(freq_separation) || freq_separation < 0
        error('Frequency separation must be nonnegative and finite.');
    end

    % Reset the selected adaptive controller
    modScheme = adaptive_controller([], method, true);    % 'BPSK', 'QPSK', '16-QAM', '64-QAM'#CHANGE
    % method = 'immediate', 'threshold', 'per_modulation',
    %          'history', 'window', 'aimd', 'score', 'q_learning', or goodput_probe'
    numBlocks = 1000;
    numSymbols = 12;

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
    % Example [10 20]: 10 ON blocks followed by 20 OFF blocks.
    % The pattern repeats continuously.
    aci_block_pattern = [10 20];

    % Mode 2: Random burst-based ACI schedule
    % [maximum duty cycle in percent, maximum consecutive ON blocks]
    % Example [40 5]: ACI is active in at most 40% of blocks,
    % with no burst longer than 5 consecutive blocks.
    randomstr_config = [40 30];

    max_power = 1;  % ACI power relative to desired-signal power

    % Amplifiers
    ampGain_Tx = 2;
    ampGain_Rx_linear = 1;
    ampGain_Rx_non_linear = 0.15;

    % Cubic nonlinearity
    a1 = 1;

    % AWGN variance is calculated per block from the requested SNR.

    %=========================================================
    %          VALIDATE ACI MODE AND GENERATE SCHEDULE
    %=========================================================

    aciSchedule = false(numBlocks, 1);

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

    %=========================================================
    %             CONSTELLATION SETUP
    %=========================================================

    switch modScheme
        case 'BPSK'
            refSyms = pskmod(0:1, 2);
        case 'QPSK'
            refSyms = pskmod(0:3, 4, pi/4);
        otherwise
            refSyms = qammod(0:M-1, M, 'UnitAveragePower', true);
    end

    % Live constellation/PSD figures are disabled for fast parameter sweeps.

    %=========================================================
    %                  SIMULATION LOOP
    %=========================================================

    fprintf('Starting %s Transmission | Frequency separation = %.2f MHz...\n', ...
        modScheme, freq_separation/1e6);
    fprintf('Block\tErrors\tStatus\n');
    fprintf('------------------------------\n');

    nackHistory = false(numBlocks, 1);
    ackHistory = false(numBlocks, 1);
    aciHistory = false(numBlocks, 1);

    current_time = 0;

    for b = 1:numBlocks

        %% ---------------------------------------------------
        % 1. TRANSMITTER
        % ---------------------------------------------------

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

        %----------------------------------------------------
        % Save the desired signal BEFORE adding ACI
        %----------------------------------------------------

        desiredPower = mean(abs(txBB).^2);

        %----------------------------------------------------
        % Adjacent-channel interference generation
        %----------------------------------------------------

        interferencePower = max_power * desiredPower;

        %----------------------------------------------------
        % RANDOM ACI MODULATION SELECTION
        %----------------------------------------------------

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

        % Save desired-only baseband waveform before ACI addition
        txBB_noACI = txBB;

        % Continuous time tracking
        t = current_time + (0:L-1)' / Fs;

        current_time = current_time + L / Fs;

        % Generate complex baseband ACI shifted in frequency
        aci_BB_raw = aci_waveform .* ...
            exp(1j * 2 * pi * freq_separation * t);

        rawAciPower = mean(abs(aci_BB_raw).^2);

        % Scale ACI to requested power
        aci_BB = sqrt(interferencePower / max(rawAciPower, eps)) ...
            .* aci_BB_raw;

        %----------------------------------------------------
        % ACI block schedule: periodic ratioed or randomstr
        %----------------------------------------------------

        switch lower(interference_type)

            case 'ratioed'

                blockPosition = mod(b - 1, aci_cycle_len);

                aciBlockOn = (blockPosition < aci_on_blocks);

            case 'randomstr'

                aciBlockOn = aciSchedule(b);

        end

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

        %% ---------------------------------------------------
        % 2. TX AMPLIFIER
        % ---------------------------------------------------

        % Desired signal after TX amplifier, BEFORE ACI
        txAmp_noACI = txBB_noACI * ampGain_Tx;

        % Actual transmitted signal after TX amplifier, WITH ACI
        txAmp = txBB * ampGain_Tx;

        %% ---------------------------------------------------
        % 3. AWGN CHANNEL
        % ---------------------------------------------------

        channelGain = 0.5;

        % SNR is measured at RX input using desired signal power only,
        % before ACI, RX gain and RX nonlinearity.
        desiredRxPower = mean(abs(channelGain * ...
            txAmp_noACI * ampGain_Rx_linear).^2);
        noiseVar = desiredRxPower / (10^(snrDb/10));

        noise = sqrt(noiseVar/2) * ...
            (randn(size(txAmp)) + 1j*randn(size(txAmp)));

        rxWaveform = channelGain * txAmp + noise;

        %% ---------------------------------------------------
        % 4. LINEAR RX AMPLIFIER
        % ---------------------------------------------------

        rxAmp = rxWaveform * ampGain_Rx_linear;

        %% ---------------------------------------------------
        % 5. NON-LINEAR RX AMPLIFIER
        % ---------------------------------------------------

        % Baseband equivalent cubic nonlinearity:
        % y = a1*x + a3*|x|^2*x

        rxNL = a1 * rxAmp + ...
            ampGain_Rx_non_linear * (abs(rxAmp).^2) .* rxAmp;

        rxWaveform = rxNL;

        % Save the waveform at the receiver input (after RX
        % amplification/nonlinearity, before matched filtering).
        rxSignalForPSD = rxWaveform;

        %% ---------------------------------------------------
        % 6. RECEIVER RRC FILTER
        % ---------------------------------------------------

        rxFiltered = upfirdn(rxWaveform, hRRC, 1, sps);

        rxSyms = rxFiltered(span + 1 : span + numSymbols);

        %% ---------------------------------------------------
        % 7. DEMODULATION
        % ---------------------------------------------------

        rxBits = decoder(rxSyms, modScheme, M);

        %% ---------------------------------------------------
        % 8. ACK / NACK
        % ---------------------------------------------------

        [isAck, numErrors] = generate_ack_nack(txBits, rxBits);

        ackHistory(b) = isAck;
        nackHistory(b) = ~isAck;

        if isAck
            status = 'ACK';
        else
            status = 'NACK';
        end

        modScheme = adaptive_controller(isAck, method, false);% 'BPSK', 'QPSK', '16-QAM', '64-QAM'#CHANGE
        % method = 'immediate', 'threshold', 'per_modulation',
        %          'history', 'window', 'aimd', 'score', 'q_learning', or goodput_probe'

        % %% ---------------------------------------------------
        % % 9. UPDATE LIVE CONSTELLATION
        % % ---------------------------------------------------
        % 
        % switch currentModScheme
        % 
        %     case 'BPSK'
        %         refSyms = pskmod(0:1, 2);
        % 
        %     case 'QPSK'
        %         refSyms = pskmod(0:3, 4, pi/4);
        % 
        %     case '16-QAM'
        %         refSyms = qammod(0:15, 16, ...
        %             'UnitAveragePower', true);
        % 
        %     case '64-QAM'
        %         refSyms = qammod(0:63, 64, ...
        %             'UnitAveragePower', true);
        % 
        % end
        % 
        % set(hIdeal, 'XData', real(refSyms), ...
        %     'YData', imag(refSyms));
        % 
        % set(hScatter, 'XData', real(rxSyms), ...
        %     'YData', imag(rxSyms));
        % 
        % maxVal = 1.5 * max(abs([real(refSyms(:)); ...
        %     imag(refSyms(:))]));
        % 
        % axis(ax, [-maxVal maxVal -maxVal maxVal]);
        % axis(ax, 'square');
        % 
        % title(ax, sprintf('%s | Block %d/%d | Errors = %d | %s', ...
        %     currentModScheme, b, numBlocks, numErrors, status));
        % 
        % %----------------------------------------------------
        % % 10. UPDATE TRANSMITTED AND RECEIVED PSD
        % %
        % % Transmitted PSD: desired signal after TX amplifier,
        % % BEFORE adding ACI.
        % %
        % % Received PSD: after channel, AWGN, RX amplifier and
        % % RX nonlinearity; includes the effects of ACI.
        % %----------------------------------------------------
        % 
        % psdWindowLength = min(2048, length(txAmp_noACI));
        % 
        % psdWindow = hann(psdWindowLength, 'periodic');
        % 
        % psdOverlap = floor(psdWindowLength / 2);
        % 
        % % Transmitter PSD: desired signal only, before ACI addition
        % [Ptx, fPSD] = pwelch(txAmp_noACI, psdWindow, psdOverlap, ...
        %     psdNFFT, Fs, 'centered', 'power');
        % 
        % % Receiver PSD: received signal after channel and RX nonlinearity
        % [Prx, ~] = pwelch(rxSignalForPSD, psdWindow, psdOverlap, ...
        %     psdNFFT, Fs, 'centered', 'power');
        % 
        % % Convert power spectrum to dBm using a 1-ohm reference
        % Ptx_dBm = 10*log10((Ptx / psdResistance) + eps) + 30;
        % Prx_dBm = 10*log10((Prx / psdResistance) + eps) + 30;
        % 
        % % Update both PSD plots
        % set(txPSDLine, 'XData', fPSD/1e6, 'YData', Ptx_dBm);
        % 
        % set(rxPSDLine, 'XData', fPSD/1e6, 'YData', Prx_dBm);
        % 
        % title(axTxPSD, sprintf('Desired TX PSD (Before ACI) | Block %d', b));
        % 
        % title(axRxPSD, sprintf('Received PSD (With ACI Effects) | Block %d', b));
        % 
        % drawnow;
        % pause(0.2);
        % 
        % if ~isvalid(hFig)
        % 
        %     fprintf('Simulation stopped by user.\n');
        %     break;
        % 
        % end

    end

    %% =====================================================
    % FINAL GOODPUT CALCULATION
    % ======================================================

    schemeNames = {'BPSK', 'QPSK', '16-QAM', '64-QAM'};

    schemeRates = [1 2 4 6];

    blockRates = zeros(completedBlocks, 1);

    for b = 1:completedBlocks

        schemeIdx = find(strcmp(schemeNames, usedModScheme{b}), 1);

        blockRates(b) = schemeRates(schemeIdx);

    end

    successfulBitsPerBlock = blockRates .* numSymbols .* ...
        double(ackHistory(1:completedBlocks));

    totalSuccessfulBits = sum(successfulBitsPerBlock);

    averageGoodput = totalSuccessfulBits / ...
        (completedBlocks * numSymbols);

    % Assign the function outputs used by the frequency-separation sweep.
    goodput = averageGoodput;
    ackRate = mean(ackHistory(1:completedBlocks));
    successfulBits = totalSuccessfulBits;

    % Modulation-history plotting is disabled for fast sweeps.

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
% Methods: immediate, threshold, per_modulation, history, window, aimd, score, goodput_probe, q_learning.
%
% ack    = 1 for ACK, 0 for NACK
% method = 'immediate', 'threshold', 'per_modulation',
%          'history', 'window', 'aimd', 'score', 'goodput_probe', or 'q_learning'
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
               'history, window, aimd, score, goodput_probe, or q_learning.']);
end


modScheme = schemes{idx};
end