function Adaptive_Comm_Sys_Baseband()
    % --- System Parameters ---
    modScheme  = adaptive_controller([], 'history', true);   
    numBlocks  = 100;
    numSymbols = 12;
    
    % Replaced RF parameters with Baseband parameters
    Rs        = 10e6;        % Symbol rate (10 MBaud)
    Bandwidth = 18e6;        % Target RF Bandwidth (15 MHz)
    sps       = 16;          % Baseband samples per symbol
    Fs        = Rs * sps;    % Baseband Sampling Frequency (160 MHz)

    % ACI parameters
    max_power       = 0.5;   % ACI/desired power ratio
    freq_separation = 15e6;  % Adjacent carrier separation = 15 MHz
    max_duty_cycle  = 0.25 ;

    % Amplifiers
    ampGain_Tx            = 2;
    ampGain_Rx_linear     = 1;
    ampGain_Rx_non_linear = 0.00;
    
    % Cubic Nonlinearity specific term
    a1 = 1;

    % AWGN variance at Room Temperature
    noiseVar = 2.055e-21;

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
    span    = 6;
    rolloff = (Bandwidth / Rs) - 1;

    hRRC = rcosdesign(rolloff, span, sps, 'sqrt');
    hRRC = hRRC / norm(hRRC); % Energy normalization from Baseband script

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
            refSyms = qammod(0:M-1, M, 'UnitAveragePower', true);
    end

    % Create figure
    hFig = figure('Name', 'Live Constellation Diagram', 'NumberTitle', 'off', 'Color', 'w');
    ax = axes(hFig);

    hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'square');

    % Axis limits
    maxVal = max(abs([real(refSyms(:)); imag(refSyms(:))])) * 1.5;
    axis(ax, [-maxVal maxVal -maxVal maxVal]);
    xlabel(ax, 'In-Phase (I)'); ylabel(ax, 'Quadrature (Q)');
    title(ax, sprintf('Live %s Constellation', modScheme));

    % Ideal constellation
    hIdeal = plot(ax, real(refSyms), imag(refSyms), 'r*', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Ideal');

    % Received constellation
    hScatter = plot(ax, NaN, NaN, 'b.', 'MarkerSize', 8, 'DisplayName', 'Received');

    legend(ax, 'Location', 'best');

    % =========================================================
    %                  SIMULATION LOOP
    % =========================================================

    fprintf('Starting %s Transmission...\n', modScheme);
    fprintf('Block\tErrors\tStatus\n');
    fprintf('------------------------------\n');
    nackHistory = false(numBlocks, 1);
    ackHistory = false(numBlocks, 1);  

    current_time = 0; % Track continuous time for complex exponent shift

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

        % RRC pulse shaping
        txWaveform = upfirdn(txSyms, hRRC, sps);
        txWaveform = txWaveform(:); % Complex baseband signal

        txBB = txWaveform;

        %----------------------------------------------------
        % Adjacent-channel interference
        %----------------------------------------------------

        desiredPower = mean(abs(txBB).^2);
        interferencePower = max_power * desiredPower;

        aci_bits = randi([0 1], blockSize, 1);
        aci_syms = encoder(aci_bits, modScheme, M);

        % RRC pulse shaping for ACI (use same filter logic as desired)
        aci_waveform = upfirdn(aci_syms, hRRC, sps);
        aci_waveform = aci_waveform(:);

        L = min(length(txBB), length(aci_waveform));
        txBB = txBB(1:L);
        aci_waveform = aci_waveform(1:L);
        
        t = current_time + (0:L-1)' / Fs;
        current_time = current_time + L / Fs;

        % Generate unscaled Complex Baseband interferer shifted to +15 MHz
        aci_BB_raw = aci_waveform .* exp(1j * 2 * pi * freq_separation * t);

        rawAciPower = mean(abs(aci_BB_raw).^2);

        % Scale to the requested average power
        aci_BB = sqrt(interferencePower / max(rawAciPower, eps)) .* aci_BB_raw;

        % Randomize the interference with a maximum duty cycle
        maxInterfered = floor(max_duty_cycle * numSymbols);
        activeSymbols = false(numSymbols, 1);
        numInterfered = randi([0, maxInterfered]);

        if numInterfered > 0
            selected = randperm(numSymbols, numInterfered);
            activeSymbols(selected) = true;
        end

        N = numel(aci_BB);
        sampleIndex = (0:N-1).';
        symbolIndex = min(floor(sampleIndex / sps) + 1, numSymbols);
        aciMask = double(activeSymbols(symbolIndex));

        aci_interference = aciMask .* aci_BB;

        % Add interference to the desired baseband signal
        txBB = txBB + aci_interference;
        
        %% ---------------------------------------------------
        % 2. TX AMPLIFIER
        % ---------------------------------------------------
        txAmp = txBB * ampGain_Tx;

        %% ---------------------------------------------------
        % 3. AWGN CHANNEL
        % ----------------------------
        % Complex Gaussian noise implementation
        noise = sqrt(noiseVar/2) * (randn(size(txAmp)) + 1j*randn(size(txAmp)));
        channelGain = 0.5;
        rxWaveform = channelGain * txAmp + noise;

        %% ---------------------------------------------------
        % 4. LINEAR RX AMPLIFIER
        % ---------------------------------------------------
        rxAmp = rxWaveform * ampGain_Rx_linear;

        %% ---------------------------------------------------
        % 5. NON-LINEAR RX AMPLIFIER
        % ---------------------------------------------------
        % Baseband Equivalent Cubic nonlinear amplifier: y = a1*x + a3*|x|^2*x
        rxNL = a1 * rxAmp + ampGain_Rx_non_linear * (abs(rxAmp).^2) .* rxAmp;
        
        rxWaveform = rxNL;

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
        
        modScheme = adaptive_controller(isAck, 'history', false); 

        fprintf('%d\t\t%d\t\t%s\n', b, numErrors, status);

        %% ---------------------------------------------------
        % 9. UPDATE LIVE CONSTELLATION
        % ---------------------------------------------------
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

        set(hIdeal, 'XData', real(refSyms), 'YData', imag(refSyms));
        set(hScatter, 'XData', real(rxSyms), 'YData', imag(rxSyms));

        maxVal = 1.5 * max(abs([real(refSyms(:)); imag(refSyms(:))]));
        axis(ax, [-maxVal maxVal -maxVal maxVal]);
        axis(ax, 'square');

        title(ax, sprintf('%s | Block %d/%d | Errors = %d | %s', currentModScheme, b, numBlocks, numErrors, status));

        drawnow;
        pause(0.2);

        if ~isvalid(hFig)
            fprintf('Simulation stopped by user.\n');
            break;
        end
    end
    
    %% FINAL GOODPUT CALCULATION
    schemeNames = {'BPSK', 'QPSK', '16-QAM', '64-QAM'};
    schemeRates = [1 2 4 6];

    blockRates = zeros(completedBlocks, 1);
    for b = 1:completedBlocks
        schemeIdx = find(strcmp(schemeNames, usedModScheme{b}), 1);
        blockRates(b) = schemeRates(schemeIdx);
    end

    successfulBitsPerBlock = blockRates .* numSymbols .* double(ackHistory(1:completedBlocks));
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
        schemeLevels(b) = find(strcmp(schemeNames, usedModScheme{b}), 1);
    end

    figure('Name', 'Adaptive Modulation History', 'NumberTitle', 'off', 'Color', 'w');

    stairs(1:completedBlocks, schemeLevels, 'b-', 'LineWidth', 1.8);
    hold on;

    nackBlocks = find(nackHistory(1:completedBlocks));
    plot(nackBlocks, schemeLevels(nackBlocks), 'rx', 'MarkerSize', 10, 'LineWidth', 2);

    grid on;
    xlim([1 max(2, completedBlocks)]);
    ylim([0.5 4.5]);

    yticks(1:4);
    yticklabels(schemeNames);

    xlabel('Block Number');
    ylabel('Modulation Scheme');
    title('Adaptive Modulation with ACK/NACK Feedback');
    legend('Modulation used', 'NACK', 'Location', 'best');
    set(gca, 'FontSize', 11);
    hold off;
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
            syms = qammod(intData, M, 'UnitAveragePower', true);
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
            rxInt = qamdemod(syms, M, 'UnitAveragePower', true);
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

function modScheme = adaptive_controller(ack, method, reset)
    persistent idx ackCount nackCount
    persistent histAck
    persistent blockCount lastSwitch
    persistent aimdScore

    persistent modHistory modObs probeCycle historyLocked
    persistent holdRemaining
    persistent scoreValue

    persistent gpHistory gpBlocksSinceProbe
    persistent gpProbeActive gpProbeBase gpProbeCount gpProbeAcks

    schemes = {'BPSK','QPSK','16-QAM','64-QAM'};
    rates = [1 2 4 6];

    if nargin < 2 || isempty(method)
        method = 'threshold';
    end
    if nargin < 3
        reset = false;
    end

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
    end

    m = 3;
    n = 1;
    mLevels = [3 3 4 5];
    nLevels = [1 1 1 1];
    historyWindow = 5;
    holdDuration = 5;
    Wwindow = 5;
    TL = 0.40;
    TH = 0.80;
    dwellTime = 2;
    alpha = 0.5;
    beta = 0.70;
    xmin = 1;
    xmax = 6;
    levels = [1 2 4 6];
    scoreA = 0.5;
    scoreB = 1.0;
    scoreMin = 0;
    scoreMax = 10;
    scoreQL = 3;
    scoreQH = 7;
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

    usedIdx = idx;
    blockCount = blockCount + 1;

    switch lower(method)
        case 'immediate'
            if ack == 1
                idx = min(idx + 1, 4);
            else
                idx = max(idx - 1, 1);
            end

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

        case 'history'
            modHistory{usedIdx} = [modHistory{usedIdx}, ack];
            if numel(modHistory{usedIdx}) > historyWindow
                modHistory{usedIdx} = modHistory{usedIdx}(end-historyWindow+1:end);
            end
            modObs(usedIdx) = min(modObs(usedIdx) + 1, historyWindow);
            if all(modObs >= historyWindow)
                historyLocked = true;
            end

            if ~historyLocked
                if usedIdx == 4
                    if ack == 0
                        probeCycle = mod(probeCycle, 3) + 1;
                        probeSchemes = [3 2 1];
                        idx = probeSchemes(probeCycle);
                    else
                        idx = 4;
                    end
                else
                    idx = 4;
                end
            else
                if usedIdx == 4 && ack == 0
                    candidates = [3 2 1];
                    goodput = zeros(1, numel(candidates));
                    for k = 1:numel(candidates)
                        s = candidates(k);
                        ackRatio = mean(modHistory{s});
                        goodput(k) = rates(s) * ackRatio;
                    end
                    [~, bestPos] = max(goodput);
                    idx = candidates(bestPos);
                    holdRemaining = holdDuration - 1;
                elseif usedIdx ~= 4 && holdRemaining > 0
                    holdRemaining = holdRemaining - 1;
                    if holdRemaining == 0
                        idx = 4;
                    else
                        idx = usedIdx;
                    end
                else
                    idx = 4;
                    holdRemaining = 0;
                end
            end

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

        case 'aimd'
            if ack == 1
                aimdScore = min(aimdScore+alpha,xmax);
            else
                aimdScore = max(beta*aimdScore,xmin);
            end
            [~,idx] = min(abs(levels-aimdScore));

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

        case 'goodput_probe'
            gpHistory{usedIdx} = [gpHistory{usedIdx}, ack];
            if numel(gpHistory{usedIdx}) > gpHistoryWindow
                gpHistory{usedIdx} = gpHistory{usedIdx}(end-gpHistoryWindow+1:end);
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
                    if probeGoodput >= baseGoodput + probeDelta && probeFailureRate <= probeEpsilon
                        idx = usedIdx;
                        fprintf(['Probe accepted: %s, goodput %.3f ' '(baseline %.3f)\n'], schemes{usedIdx}, probeGoodput, baseGoodput);
                    else
                        idx = gpProbeBase;
                        fprintf(['Probe rejected: %s, goodput %.3f ' '(baseline %.3f), failure rate %.1f%%\n'], schemes{usedIdx}, probeGoodput, baseGoodput, 100 * probeFailureRate);
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
                if gpBlocksSinceProbe >= probeInterval && usedIdx < 4 && ~isempty(gpHistory{usedIdx})
                    gpProbeBase = usedIdx;
                    gpProbeActive = true;
                    gpProbeCount = 0;
                    gpProbeAcks = 0;
                    idx = usedIdx + 1;
                end
            end
        otherwise
            error(['Unknown method. Use immediate, threshold, per_modulation, ' 'history, window, aimd, score, or goodput_probe.']);
    end

    modScheme = schemes{idx};
end