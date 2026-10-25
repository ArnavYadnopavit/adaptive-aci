function Generic_Comm_Sim_With_Plot()

    % --- System Parameters ---
    modScheme  = adaptive_controller([], 'goodput_probe', true);   % 'BPSK', 'QPSK', '16-QAM', '64-QAM'#CHANGE
    % method = 'immediate', 'threshold', 'per_modulation',
    %          'history', 'window', 'aimd', 'score', or 'goodput_probe'
    numBlocks  = 400;
    numSymbols = 12;
    fc = 2.4e9;   % Carrier frequency
    Fs = 50e9;    % Sampling frequency
    Rs = 1e6;     % Symbol rate

    %ACI parameters
    interference_type = 'ratioed';
    ratio = [10 20]; %consecutive blocks without aci to aci ratio
    max_power = 0.05 ;
    freq_seperation = 1.25e3/2; 
    adjacentCenter = fc + freq_seperation;

    maxDutyCycle = 0.40;       %maximum fraction of blocks with ACI
    maxBurstLength = 5;        %maximum consecutive ACI blocks

    % ---- ACI schedule (only used when interference_type = 'random') ----
    aciSchedule = false(numBlocks, 1);
    if strcmpi(interference_type, 'random')
        aciSchedule = random_aci_schedule(numBlocks, 4);
    end
    % Generate ACI schedule for random and randomstr modes
    aciSchedule = false(numBlocks, 1);

    if strcmpi(interference_type, 'random')
        aciSchedule = random_aci_schedule(numBlocks, 4);

    elseif strcmpi(interference_type, 'randomstr')

        maxInterferenceBlocks = floor(maxDutyCycle * numBlocks);

        bSch = 1;
        interferenceBlocksUsed = 0;

        while bSch <= numBlocks && ...
                interferenceBlocksUsed < maxInterferenceBlocks

            if rand < 0.5
                remainingBlocks = numBlocks - bSch + 1;
                remainingBudget = maxInterferenceBlocks ...
                    - interferenceBlocksUsed;

                burstLength = randi([1, min([remainingBlocks, ...
                    remainingBudget, maxBurstLength])]);

                aciSchedule(bSch:bSch+burstLength-1) = true;

                interferenceBlocksUsed = interferenceBlocksUsed ...
                    + burstLength;

                cleanGap = randi([1, ...
                    max(1, numBlocks-bSch-burstLength+1)]);

                bSch = bSch + burstLength + cleanGap;
            else
                bSch = bSch + 1;
            end
        end
    end


    % Amplifiers
    ampGain_Tx          = 2;
    %% 
    ampGain_Rx_linear   = 1;
    ampGain_Rx_non_linear = 0.09;

    % AWGN variance
    noiseVar = 1e-3;

   
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
    sps     = 50000;
    span    = 6;
    rolloff = 0.25;

    hRRC = rcosdesign(rolloff, span, sps, 'sqrt');

    %-----ACI parameters --------

    sps_aci     = 50000;
    span_aci    = 6;
    rolloff_aci = 0.25;

    hRRC_aci = rcosdesign(rolloff_aci, span, sps, 'sqrt');





    %=========================================================
    %             CONSTELLATION SETUP
    % =========================================================

    % Generate ideal constellation points
    switch modScheme

        case 'BPSK'
            refSyms = pskmod(0:1, 2);

        case 'QPSK'
            refSyms = pskmod(0:3, 4, pi/4);

        otherwise
            refSyms = qammod(0:M-1, M, ...
                'UnitAveragePower', true);
    end


    % Create figure
    hFig = figure( ...
        'Name', 'Live Constellation Diagram', ...
        'NumberTitle', 'off', ...
        'Color', 'w');

    ax = axes(hFig);

    hold(ax, 'on');
    grid(ax, 'on');
    axis(ax, 'square');

    % Axis limits
    maxVal = max(abs([real(refSyms(:)); ...
        imag(refSyms(:))])) * 1.5;

    axis(ax, [-maxVal maxVal -maxVal maxVal]);

    xlabel(ax, 'In-Phase (I)');
    ylabel(ax, 'Quadrature (Q)');

    title(ax, sprintf( ...
        'Live %s Constellation', modScheme));


    % Ideal constellation
    hIdeal = plot(ax, ...
        real(refSyms), ...
        imag(refSyms), ...
        'r*', ...
        'MarkerSize', 10, ...
        'LineWidth', 1.5, ...
        'DisplayName', 'Ideal');

    % Received constellation
    hScatter = plot(ax, ...
        NaN, NaN, ...
        'b.', ...
        'MarkerSize', 8, ...
        'DisplayName', 'Received');


    legend(ax, 'Location', 'best');



    % =========================================================
    %                  SIMULATION LOOP
    % =========================================================

    fprintf('Starting %s Transmission...\n', modScheme);

    fprintf('Block\tErrors\tStatus\n');
    fprintf('------------------------------\n');
    nackHistory = false(numBlocks, 1);
    ackHistory = false(numBlocks, 1);  % ACK/NACK outcome for each block
    aciHistory = zeros(1, numBlocks);

    for b = 1:numBlocks

        %% ---------------------------------------------------
        % 1. TRANSMITTER
        % ---------------------------------------------------
        % Record the modulation used for this block
        currentModScheme = modScheme;
        usedModScheme{b} = currentModScheme;
        completedBlocks = b;

        % Update modulation parameters for this block
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

        % RRC pulse shaping
        txWaveform = upfirdn( ...
            txSyms, ...
            hRRC, ...
            sps);
        % Plot the pulse-shaped TX waveform for the first block
        n = (0:numel(txWaveform)-1).';
        t = n / Fs;
        
        txRF = sqrt(2) * real( ...
            txWaveform(:) .* exp(1j*2*pi*fc*t));

       


        %----------------------------------------------------
        % Adjacent-channel interference
        %----------------------------------------------------


        % Desired RF power BEFORE adding interference
        desiredPower = mean(abs(txRF(:)).^2);

        % Here max_power is an interference-to-desired POWER RATIO
        % max_power = 1 means equal unmasked average powers
        interferencePower = max_power * desiredPower;

        aci_bits = randi([0 1], blockSize, 1);
        aci_syms = encoder(aci_bits, modScheme, M);


        % RRC pulse shaping
        aci_waveform = upfirdn( ...
            aci_syms, ...
            hRRC_aci, ...
            sps);

        % Generate unscaled RF interferer
        aci_RF_raw = sqrt(2) * real( ...
            aci_waveform(:) .* exp(1j*2*pi*adjacentCenter*t));

        % Measure its existing power
        rawAciPower = mean(abs(aci_RF_raw).^2);

        % Scale to the requested average power
        aci_RF = sqrt(interferencePower / max(rawAciPower, eps)) ...
            .* aci_RF_raw;
        switch lower(interference_type)

            
case 'ratioed'

    cleanBlocks = ratio(1);
    aciBlocks   = ratio(2);

    if cleanBlocks < 0 || aciBlocks < 1 || ...
            cleanBlocks ~= floor(cleanBlocks) || ...
            aciBlocks ~= floor(aciBlocks)
        error('ratio must be [nonnegative cleanBlocks, positive aciBlocks].');
    end

    cycleLength = cleanBlocks + aciBlocks;
    positionInCycle = mod(b-1, cycleLength);

    % First cleanBlocks have no ACI; next aciBlocks have ACI
    aciEnabled = positionInCycle >= cleanBlocks;

    if aciEnabled
        aci_interference = aci_RF;
    else
        aci_interference = zeros(size(aci_RF));
    end
            case 'random'
                aciEnabled = aciSchedule(b);

                if aciEnabled
                    aci_interference = aci_RF;
                else
                    aci_interference = zeros(size(aci_RF));
                end
            case 'randomstr'
                aciEnabled = aciSchedule(b);

                if aciEnabled
                    aci_interference = aci_RF;
                else
                    aci_interference = zeros(size(aci_RF));
                end
            otherwise
                error('Unknown interference_type: %s', interference_type);

        end

        % Add interference to the desired RF signal
        txRF = txRF + aci_interference;

        % Display the number of interfered symbols
        fprintf('Block %d: ACI = %d\n', b, aciEnabled);
        
        %% ---------------------------------------------------
        % 2. TX AMPLIFIER
        % ---------------------------------------------------

        txAmp = txWaveform * ampGain_Tx;
        txAmpRF = txRF*ampGain_Tx ;

        %% ---------------------------------------------------
        % 3. AWGN CHANNEL
        % ----------------------------
        rxWaveform = awgn_channel( ...
            txAmp, ...
            noiseVar);
       rxWaveformRF = awgn_channel_rf(txAmpRF,noiseVar) ;


        %% ---------------------------------------------------
        % 4. LINEAR RX AMPLIFIER
        % ---------------------------------------------------

        rxAmp = rxWaveform * ampGain_Rx_linear;
        rxAmpRF = rxWaveformRF * ampGain_Rx_linear ;

        %% ---------------------------------------------------
        % 5. NON-LINEAR RX AMPLIFIER
        % ---------------------------------------------------

        rxNL = rxAmp + ...
               ampGain_Rx_non_linear * (rxAmp.^3);


        rxWaveform = rxNL;
        rxWaveformRF = rxAmpRF + 0 ;

        %% ---------------------------------------------------
        % Demodulating after ACI
        % ---------------------------------------------------
        t = (0:numel(rxWaveformRF)-1).' / Fs;

        rxBasebandRF = sqrt(2) * rxWaveformRF(:) .* ...
            exp(-1j*2*pi*fc*t);

        % RRC matched filtering also removes the high-frequency mixing image
        rxFilteredRF = upfirdn(rxBasebandRF, hRRC, 1, sps);

        % Remove combined transmit/receive RRC delay
        rxSymsRF = rxFilteredRF(span + 1 : span + numSymbols);





        %% ---------------------------------------------------
        % 6. RECEIVER RRC FILTER
        % ---------------------------------------------------

        rxFiltered = upfirdn( ...
            rxWaveform, ...
            hRRC, ...
            1, ...
            sps);
       

        % Remove RRC filter delay
        rxSyms = rxFiltered( ...
            span + 1 : ...
            span + numSymbols);
       


        % Undo TX gain
        


        %% ---------------------------------------------------
        % 7. DEMODULATION
        % ---------------------------------------------------

        rxBits = decoder( ...
            rxSymsRF, ...
            modScheme, ...
            M);


        %% ---------------------------------------------------
        % 8. ACK / NACK
        % ---------------------------------------------------

        [isAck, numErrors] = ...
            generate_ack_nack(txBits, rxBits);
        ackHistory(b) = isAck;
        nackHistory(b) = ~isAck;

        if isAck
            status = 'ACK';
        else
            status = 'NACK';
        end
        modScheme = adaptive_controller(isAck, 'goodput_probe', false); %#CHANGE


        fprintf('%d\t\t%d\t\t%s\n', ...
            b, numErrors, status);



        %% ---------------------------------------------------
        % 9. UPDATE LIVE CONSTELLATION
        % ---------------------------------------------------

        % Use the modulation scheme used for THIS block
        %{
        switch currentModScheme
            case 'BPSK'
                refSyms = pskmod(0:1, 2);

            case 'QPSK'
                refSyms = pskmod(0:3, 4, pi/4);

            case '16-QAM'
                refSyms = qammod(0:15, 16, 'UnitAveragePower', true);

            case '64-QAM'
                refSyms = qammod(0:63, 64, 'UnitAveragePower', true);
        end

        % Update ideal constellation points
        % Ideal constellation — SAVE the plot handle
        
        % Update received symbols from the same receiver path used for decoding
        set(hIdeal, ...
    'XData', real(refSyms), ...
    'YData', imag(refSyms));
        set(hScatter, ...
            'XData', real(rxSymsRF), ...
            'YData', imag(rxSymsRF));

        % Rescale axes
        maxVal = 1.5 * max(abs([real(refSyms(:)); imag(refSyms(:))]));
        axis(ax, [-maxVal maxVal -maxVal maxVal]);
        axis(ax, 'square');

        title(ax, sprintf( ...
            '%s | Block %d/%d | Errors = %d | %s', ...
            currentModScheme, b, numBlocks, numErrors, status));

        drawnow;
        pause(0.2);

        if ~isvalid(hFig)
            fprintf('Simulation stopped by user.\n');
            break;
        end
        %}
        aciHistory(b) = aciEnabled;

    end
    %% FINAL GOODPUT CALCULATION

    % Payload rates in bits/symbol: BPSK, QPSK, 16-QAM, 64-QAM.
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
    averageGoodput = totalSuccessfulBits / (completedBlocks * numSymbols);
    ackRate = mean(ackHistory(1:completedBlocks));

    fprintf('\n========== FINAL SIMULATION RESULTS ==========\n');
    fprintf('Completed blocks        : %d\n', completedBlocks);
    fprintf('ACK rate                : %.2f%%\n', 100 * ackRate);
    fprintf('Total successful payload: %.0f bits\n', totalSuccessfulBits);
    fprintf('Average goodput         : %.4f bits/symbol\n', averageGoodput);
    fprintf('==============================================\n\n');
    %% MODULATION SCHEME VS BLOCK NUMBER

    schemeNames = {'BPSK', 'QPSK', '16-QAM', '64-QAM'};
    schemeLevels = zeros(completedBlocks, 1);

    for b = 1:completedBlocks
        schemeLevels(b) = find(strcmp( ...
            schemeNames, usedModScheme{b}), 1);
    end

    figure('Name', 'Adaptive Modulation History', ...
        'NumberTitle', 'off', 'Color', 'w');

    stairs(1:completedBlocks, schemeLevels, ...
        'b-', 'LineWidth', 1.8);
    hold on;

    % Mark NACKs with red crosses at the modulation level used
    nackBlocks = find(nackHistory(1:completedBlocks));

    plot(nackBlocks, schemeLevels(nackBlocks), ...
        'rx', 'MarkerSize', 10, 'LineWidth', 2);

    grid on;
    xlim([1 max(2, completedBlocks)]);
    ylim([0.5 4.5]);

    yticks(1:4);
    yticklabels(schemeNames);

    xlabel('Block Number');
    ylabel('Modulation Scheme');
    % Right axis: ACI status
yyaxis right
stairs(1:completedBlocks, aciHistory(1:completedBlocks), ...
       'k-', 'LineWidth', 1.8);

ylim([-0.1 1.1]);
yticks([0 1]);
yticklabels({'ACI OFF', 'ACI ON'});
ylabel('ACI Status');

grid on;
xlim([1 max(2, completedBlocks)]);

xlabel('Block Number');
title('Adaptive Modulation with ACK/NACK and ACI');

legend('Modulation used', 'NACK', 'ACI Status', ...
       'Location', 'best');

set(gca, 'FontSize', 11);
hold off;

    legend('Modulation used', 'NACK', ...
        'Location', 'best');

    set(gca, 'FontSize', 11);
    hold off;

end


%% =========================================================
%                    ENCODER
% ==========================================================

function syms = encoder(bits, modScheme, M)

    k = log2(M);

    bitMatrix = reshape(bits, k, [])';

    intData = bi2de( ...
        bitMatrix, ...
        'left-msb');


    switch modScheme

        case 'BPSK'
            syms = pskmod(intData, 2);

        case 'QPSK'
            syms = pskmod( ...
                intData, ...
                4, ...
                pi/4);

        otherwise
            syms = qammod( ...
                intData, ...
                M, ...
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
            rxInt = pskdemod( ...
                syms, ...
                4, ...
                pi/4);

        otherwise
            rxInt = qamdemod( ...
                syms, ...
                M, ...
                'UnitAveragePower', true);

    end


    k = log2(M);

    bitMatrix = de2bi( ...
        rxInt, ...
        k, ...
        'left-msb');

    bits = reshape( ...
        bitMatrix', ...
        [], ...
        1);

end

%% =========================================================
%                RANDOM ACI SCHEDULE
% ==========================================================

function aciSchedule = random_aci_schedule(numBlocks, numParts)
% Returns a logical vector (numBlocks x 1), true where ACI is active.
% The blocks are split into numParts parts, and each part is split into
% two halves. In each part, ACI starts at a uniform random block in the
% first half and stops at a uniform random block in the second half.
% Each part is drawn independently.

if nargin < 2
    numParts = 4;
end

assert(mod(numBlocks, numParts*2) == 0, ...
    'numBlocks must be divisible by %d for the random scheme.', ...
    numParts*2);

aciSchedule = false(numBlocks, 1);

partLen = numBlocks / numParts;   % e.g. 40
halfLen = partLen / 2;            % e.g. 20

for p = 1:numParts
    offset = (p-1) * partLen;

    % Start: uniform random block in the FIRST half of this part
    startBlk = offset + randi(halfLen);

    % Stop: uniform random block in the SECOND half of this part
    stopBlk  = offset + halfLen + randi(halfLen);

    % ACI is on from startBlk through stopBlk (inclusive)
    aciSchedule(startBlk:stopBlk) = true;
end
end



%% =========================================================
%                    AWGN CHANNEL
% ==========================================================

function rxWaveform = awgn_channel( ...
            txWaveform, noiseVar)

    % Complex Gaussian noise
    noise = sqrt(noiseVar/2) .* ...
        (randn(size(txWaveform)) + ...
         1j * randn(size(txWaveform)));
    h = 0.5;
    rxWaveform = h*txWaveform + noise;

end

function rxSignal = awgn_channel_rf(txSignal, noiseVar)

% Generate real Gaussian noise
noise = sqrt(noiseVar) * randn(size(txSignal));

% Channel attenuation
h = 0.5;

% Received real signal
rxSignal = h * txSignal + noise;

end

%% =========================================================
%                    ACK / NACK
% ==========================================================

function [isAck, errors] = ...
            generate_ack_nack(txBits, rxBits)

    errors = sum(txBits ~= rxBits);

    if errors == 0
        isAck = true;
    else
        isAck = false;
    end

end


function modScheme = adaptive_controller(ack, method, reset)
% Adaptive modulation using ACK/NACK feedback.
% Methods: immediate, threshold, per_modulation, history, window, aimd, score, goodput_probe.
%
% ack    = 1 for ACK, 0 for NACK
% method = 'immediate', 'threshold', 'per_modulation',
%          'history', 'window', or 'aimd'
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

schemes = {'BPSK','QPSK','16-QAM','64-QAM'};

% Relative useful payload rates (bits/symbol)
rates = [1 2 4 6];

if nargin < 2 || isempty(method)
    method = 'threshold';
end
if nargin < 3
    reset = false;
end

% Initialize/reset
if isempty(idx) || reset
    idx = 4;
    if strcmpi(method, 'goodput_probe')
        idx = 1;  % Start probing from BPSK
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
    otherwise
        error(['Unknown method. Use immediate, threshold, per_modulation, ' ...
               'history, window, aimd, score, or goodput_probe.']);
end


modScheme = schemes{idx};

end