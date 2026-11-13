function experiment5_duty_cycle_sweep()

    %=========================================================
    % SYSTEM PARAMETERS
    %=========================================================

    numBlocks = 1000;
    numSymbols = 12;

    Rs = 10e6;
    Bandwidth = 18e6;
    sps = 16;
    Fs = Rs * sps;

    %=========================================================
    % ACI PARAMETERS
    %=========================================================

    interference_type = 'ratioed';

    % Sweep ACI ON blocks from 0 to 40.
    % ACI OFF blocks = 40 - ACI ON blocks.
    aci_on_values = 0:40;
    aci_period = 40;

    % Fixed ACI power relative to desired signal power
    max_power = 1;

    % Adaptive modulation algorithm
    adaptive_method = 'history';

    freq_separation = 10e6;

    % Amplifiers
    ampGain_Tx = 2;
    ampGain_Rx_linear = 1;
    ampGain_Rx_non_linear = 0.15;

    % Cubic nonlinearity
    a1 = 1;

    % AWGN variance
    noiseVar = 2.055e-21;

    channelGain = 0.5;

    %=========================================================
    % RRC FILTER
    %=========================================================

    span = 6;
    rolloff = (Bandwidth / Rs) - 1;
    hRRC = rcosdesign(rolloff, span, sps, 'sqrt');
    hRRC = hRRC / norm(hRRC);

    %=========================================================
    % MODULATION AND GOODPUT SETTINGS
    %=========================================================

    schemeNames = {'BPSK', 'QPSK', '16-QAM', '64-QAM'};
    schemeRates = [1 2 4 6];

    dutyCyclePercent = 100 * aci_on_values / aci_period;
    aci_off_values = aci_period - aci_on_values;

    averageGoodput = zeros(size(aci_on_values));
    ackRate = zeros(size(aci_on_values));

    % Same random realization at every duty-cycle point
    baseSeed = 100;

    %=========================================================
    % DUTY-CYCLE SWEEP
    %=========================================================

    fprintf('ACI type = %s\n', interference_type);
    fprintf('Adaptive method = %s\n', adaptive_method);
    fprintf('Fixed max_power = %.2f\n', max_power);

    for p = 1:length(aci_on_values)

        aci_on_blocks = aci_on_values(p);
        aci_off_blocks = aci_off_values(p);
        aci_cycle_len = aci_period;

        %-----------------------------------------------------
        % GENERATE ACI SCHEDULE FOR THIS DUTY CYCLE
        %-----------------------------------------------------

        rng(10);
        aciSchedule = false(numBlocks, 1);

        for bSch = 1:numBlocks
            blockPosition = mod(bSch - 1, aci_cycle_len);
            aciSchedule(bSch) = blockPosition < aci_on_blocks;
        end

        actualDutyCycle = 100 * mean(aciSchedule);

        % Reset adaptive controller for every duty-cycle point
        modScheme = adaptive_controller([], adaptive_method, true);

        ackHistory = false(numBlocks, 1);
        usedModScheme = cell(numBlocks, 1);

        % Repeat the same random sequence for a fair comparison
        rng(baseSeed);

        current_time = 0;

        fprintf('\nRunning ACI ON/OFF = %d/%d, duty cycle = %.2f%% (%d/%d)\n', ...
            aci_on_blocks, aci_off_blocks, actualDutyCycle, ...
            p, length(aci_on_values));

        %-----------------------------------------------------
        % SIMULATION OVER ALL BLOCKS
        %-----------------------------------------------------

        for b = 1:numBlocks

            %---------------- TRANSMITTER --------------------

            currentModScheme = modScheme;
            usedModScheme{b} = currentModScheme;

            switch currentModScheme
                case 'BPSK'
                    M = 2;
                case 'QPSK'
                    M = 4;
                case '16-QAM'
                    M = 16;
                case '64-QAM'
                    M = 64;
                otherwise
                    error('Unsupported modulation scheme.');
            end

            k = log2(M);
            blockSize = numSymbols * k;
            txBits = randi([0 1], blockSize, 1);
            txSyms = encoder(txBits, currentModScheme, M);

            txWaveform = upfirdn(txSyms, hRRC, sps);
            txBB = txWaveform(:);

            % Desired signal power before adding ACI
            desiredPower = mean(abs(txBB).^2);

            % Fixed ACI power relative to desired signal power.
            % max_power = 1 means ACI power equals desired signal power
            % during ACI-active blocks.
            interferencePower = max_power * desiredPower;

            %---------------- ACI GENERATION -----------------

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

            aci_waveform = upfirdn(aci_syms, hRRC, sps);
            aci_waveform = aci_waveform(:);

            L = min(length(txBB), length(aci_waveform));
            txBB = txBB(1:L);
            aci_waveform = aci_waveform(1:L);

            t = current_time + (0:L-1)' / Fs;
            current_time = current_time + L / Fs;

            aci_BB_raw = aci_waveform .* ...
                exp(1j * 2 * pi * freq_separation * t);

            rawAciPower = mean(abs(aci_BB_raw).^2);

            aci_BB = sqrt(interferencePower / max(rawAciPower, eps)) ...
                .* aci_BB_raw;

            % Apply the ACI schedule for this duty-cycle point
            aciBlockOn = aciSchedule(b);

            if aciBlockOn
                aciMask = ones(L, 1);
            else
                aciMask = zeros(L, 1);
            end

            aci_interference = aciMask .* aci_BB;
            txBB = txBB + aci_interference;

            %---------------- TX AMPLIFIER -------------------

            txAmp = txBB * ampGain_Tx;

            %---------------- AWGN CHANNEL ------------------

            noise = sqrt(noiseVar/2) * ...
                (randn(size(txAmp)) + 1j * randn(size(txAmp)));

            rxWaveform = channelGain * txAmp + noise;

            %---------------- RX AMPLIFIERS -----------------

            rxAmp = rxWaveform * ampGain_Rx_linear;

            rxNL = a1 * rxAmp + ...
                ampGain_Rx_non_linear * (abs(rxAmp).^2) .* rxAmp;

            %---------------- RECEIVER FILTER ----------------

            rxFiltered = upfirdn(rxNL, hRRC, 1, sps);
            rxSyms = rxFiltered(span + 1 : span + numSymbols);

            %---------------- DEMODULATION ------------------

            rxBits = decoder(rxSyms, currentModScheme, M);

            %---------------- ACK / NACK --------------------

            [isAck, numErrors] = generate_ack_nack(txBits, rxBits);
            ackHistory(b) = isAck;

            % Update adaptive modulation for the next block
            modScheme = adaptive_controller(isAck, adaptive_method, false);

        end

        %-----------------------------------------------------
        % GOODPUT CALCULATION
        %-----------------------------------------------------

        blockRates = zeros(numBlocks, 1);

        for b = 1:numBlocks
            schemeIdx = find(strcmp(schemeNames, usedModScheme{b}), 1);
            blockRates(b) = schemeRates(schemeIdx);
        end

        successfulBitsPerBlock = blockRates .* numSymbols .* ...
            double(ackHistory);

        totalSuccessfulBits = sum(successfulBitsPerBlock);

        averageGoodput(p) = totalSuccessfulBits / ...
            (numBlocks * numSymbols);

        ackRate(p) = mean(ackHistory);

        fprintf('ACK rate        = %.2f%%\n', 100 * ackRate(p));
        fprintf('Average goodput = %.4f bits/symbol\n', ...
            averageGoodput(p));

    end

    %=========================================================
    % GOODPUT vs DUTY-CYCLE PLOT
    %=========================================================

    figure('Name', sprintf('Duty Cycle vs Goodput | Adaptive Scheme: %s', ...
        adaptive_method), ...
        'NumberTitle', 'off', 'Color', 'w');

    plot(dutyCyclePercent, averageGoodput, 'bo-', ...
        'LineWidth', 1.8, 'MarkerSize', 4, ...
        'MarkerFaceColor', 'b');

    grid on;
    xlabel('ACI Duty Cycle (%)');
    ylabel('Average Goodput (bits/symbol)');
    title(sprintf('Average Goodput vs ACI Duty Cycle (Scheme: %s, max power = %.1f)', ...
        adaptive_method, max_power));

    xlim([0 100]);
    xticks(0:10:100);

    %=========================================================
    % PRINT RESULTS
    %=========================================================

    fprintf('\n========== DUTY-CYCLE SWEEP RESULTS ==========\n');
    fprintf('ACI ON\tACI OFF\tDuty Cycle (%%)\tAverage Goodput\tACK Rate\n');
    fprintf('----------------------------------------------------------------\n');

    for p = 1:length(aci_on_values)
        fprintf('%d\t%d\t%.2f\t\t%.4f\t\t%.2f%%\n', ...
            aci_on_values(p), aci_off_values(p), dutyCyclePercent(p), ...
            averageGoodput(p), 100 * ackRate(p));
    end

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
% Methods: immediate, threshold, per_modulation, history,
%          window, aimd, score, goodput_probe, q_learning.
%
% ack    = 1 for ACK, 0 for NACK
% method = adaptive modulation algorithm
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

% Q-Learning state
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

% Q-learning parameters
qGamma    = 0.60;
qAlphaMin = 0.10;
qEps0     = 0.15;
qEpsMin   = 0.02;
qEpsDecay = 0.99;
qMaxSince = 3;

% Initialize/reset
if isempty(idx) || reset

    idx = 4;

    if strcmpi(method, 'goodput_probe')
        idx = 1;
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
        modHistory{k} = zeros(1,0);
    end

    bestModIdx = 3;

    probeCycle = 0;
    historyLocked = false;
    holdRemaining = 0;

    scoreValue = 5;

    gpHistory = cell(1,4);

    for k = 1:4
        gpHistory{k} = [];
    end

    gpBlocksSinceProbe = 0;
    gpProbeActive = false;
    gpProbeBase = 1;
    gpProbeCount = 0;
    gpProbeAcks = 0;

    % Q-learning initialization
    qTable = zeros(4*(qMaxSince+1), 4);

    for a = 1:4
        qTable(:,a) = rates(a)/(1-qGamma);
    end

    qCount = zeros(size(qTable));
    qSince = qMaxSince;
    qPrevState = qMaxSince*4 + 4;

    if strcmpi(method, 'q_learning')
        idx = 4;
    end

end

% General parameters
m = 3;
n = 1;

% Per-modulation thresholds
mLevels = [3 3 4 5];
nLevels = [1 1 1 1];

% History-based parameters
historyWindow = 5;
holdDuration = 5;

% Sliding-window hysteresis
Wwindow = 5;
TL = 0.40;
TH = 0.80;
dwellTime = 2;

% AIMD parameters
alpha = 0.5;
beta = 0.70;
xmin = 1;
xmax = 6;
levels = [1 2 4 6];

% Score controller
scoreA = 0.5;
scoreB = 1.0;
scoreMin = 0;
scoreMax = 10;
scoreQL = 3;
scoreQH = 7;

% Goodput-aware probing
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

% Save the modulation used for the current block
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

        % Record ACK/NACK for the modulation used in this block
        modHistory{usedIdx} = [modHistory{usedIdx}, ack];

        % Retain only the most recent historyWindow observations
        if numel(modHistory{usedIdx}) > historyWindow
            modHistory{usedIdx} = ...
                modHistory{usedIdx}(end-historyWindow+1:end);
        end

        modObs(usedIdx) = min(modObs(usedIdx) + 1, historyWindow);

        % Lock once every modulation has enough observations
        if ~historyLocked && all(modObs >= historyWindow)

            candidates = [3 2 1];
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

        % Phase 1: Cyclic exploration
        if ~historyLocked

            if usedIdx == 4

                if ack == 0

                    probeCycle = mod(probeCycle, 3) + 1;

                    % Cycle: 16-QAM -> QPSK -> BPSK
                    probeSchemes = [3 2 1];
                    idx = probeSchemes(probeCycle);

                else
                    idx = 4;
                end

            else
                idx = 4;
            end

        % Phase 2: Locked operation
        else

            if usedIdx == 4 && ack == 0

                idx = bestModIdx;
                holdRemaining = holdDuration - 1;

            elseif usedIdx ~= 4 && holdRemaining > 0

                holdRemaining = holdRemaining - 1;

                if holdRemaining == 0
                    idx = 4;
                else
                    idx = bestModIdx;
                end

            else

                idx = 4;
                holdRemaining = 0;

            end
        end

    %--------------------------------------------------------------
    % Algorithm 5: Sliding-window ACK ratio with hysteresis
    %--------------------------------------------------------------

    case 'window'

        histAck = [histAck ack];

        if numel(histAck) > Wwindow
            histAck = histAck(end-Wwindow+1:end);
        end

        if numel(histAck) >= Wwindow

            pACK = mean(histAck);

            canSwitch = (blockCount-lastSwitch) >= dwellTime;

            newIdx = idx;

            if canSwitch && pACK >= TH
                newIdx = min(idx+1, 4);

            elseif canSwitch && pACK <= TL
                newIdx = max(idx-1, 1);
            end

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

        gpHistory{usedIdx} = [gpHistory{usedIdx}, ack];

        if numel(gpHistory{usedIdx}) > gpHistoryWindow
            gpHistory{usedIdx} = ...
                gpHistory{usedIdx}(end-gpHistoryWindow+1:end);
        end

        if gpProbeActive && usedIdx == gpProbeBase + 1

            gpProbeCount = gpProbeCount + 1;
            gpProbeAcks = gpProbeAcks + ack;

            if gpProbeCount >= probeDuration

                probeAckRatio = gpProbeAcks / gpProbeCount;
                probeFailureRate = 1 - probeAckRatio;
                probeGoodput = rates(usedIdx) * probeAckRatio;

                baseHistory = gpHistory{gpProbeBase};

                if isempty(baseHistory)
                    baseAckRatio = 0;
                else
                    baseAckRatio = mean(baseHistory);
                end

                baseGoodput = rates(gpProbeBase) * baseAckRatio;

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
                idx = usedIdx;
            end

        elseif gpProbeActive

            idx = gpProbeBase;
            gpProbeActive = false;
            gpBlocksSinceProbe = 0;

        else

            gpBlocksSinceProbe = gpBlocksSinceProbe + 1;
            idx = usedIdx;

            if gpBlocksSinceProbe >= probeInterval && usedIdx < 4 && ...
                    ~isempty(gpHistory{usedIdx})

                gpProbeBase = usedIdx;
                gpProbeActive = true;
                gpProbeCount = 0;
                gpProbeAcks = 0;

                idx = usedIdx + 1;
            end

        end

    %--------------------------------------------------------------
    % Algorithm 9: Model warm-started belief-state Q-Learning
    %--------------------------------------------------------------

    case 'q_learning'

        if ack == 0
            qSince = 0;
        else
            qSince = min(qSince + 1, qMaxSince);
        end

        nextState = qSince*4 + usedIdx;

        reward = rates(usedIdx) * ack;

        qCount(qPrevState, usedIdx) = qCount(qPrevState, usedIdx) + 1;

        alphaQ = max(qAlphaMin, ...
            1/(1 + qCount(qPrevState, usedIdx)));

        target = reward + qGamma * max(qTable(nextState,:));

        qTable(qPrevState, usedIdx) = qTable(qPrevState, usedIdx) + ...
            alphaQ * (target - qTable(qPrevState, usedIdx));

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
