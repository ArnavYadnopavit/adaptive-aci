function Adaptive_baseband_updated_aci_with_psd()

    % --- System Parameters ---
    modScheme = adaptive_controller([], 'immediate', true);
    numBlocks = 100;
    numSymbols = 12;

    % Baseband parameters
    Rs = 10e6;
    Bandwidth = 18e6;
    sps = 16;
    Fs = Rs * sps;

    % ACI parameters
    % Pattern = [ACI_ON_blocks, ACI_OFF_blocks]
    % Example [10 20]: 10 ON blocks followed by 20 OFF blocks.
    % The pattern repeats continuously.
    aci_block_pattern = [10 20];

    max_power = 0.75;
    freq_separation = 10e6;

    % Amplifiers
    ampGain_Tx = 2;
    ampGain_Rx_linear = 1;
    ampGain_Rx_non_linear = 0.15;

    % Cubic nonlinearity
    a1 = 1;

    % AWGN variance
    noiseVar = 2.055e-21;

    % Validate ACI block pattern
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

    % Create one live figure containing:
    % 1. Received constellation
    % 2. Transmitted-signal PSD
    % 3. Received-signal PSD
    hFig = figure('Name', 'Live Constellation and PSD', ...
        'NumberTitle', 'off', 'Color', 'w');

    tiledlayout(hFig, 1, 3, 'TileSpacing', 'compact', ...
        'Padding', 'compact');

    %---------------- CONSTELLATION AXIS ----------------
    ax = nexttile;
    hold(ax, 'on');
    grid(ax, 'on');
    axis(ax, 'square');

    maxVal = max(abs([real(refSyms(:)); imag(refSyms(:))])) * 1.5;
    axis(ax, [-maxVal maxVal -maxVal maxVal]);

    xlabel(ax, 'In-Phase (I)');
    ylabel(ax, 'Quadrature (Q)');
    title(ax, sprintf('Live %s Constellation', modScheme));

    hIdeal = plot(ax, real(refSyms), imag(refSyms), ...
        'r*', 'MarkerSize', 10, 'LineWidth', 1.5, ...
        'DisplayName', 'Ideal');

    hScatter = plot(ax, NaN, NaN, ...
        'b.', 'MarkerSize', 8, 'DisplayName', 'Received');

    legend(ax, 'Location', 'best');

    %---------------- TRANSMITTED PSD AXIS ----------------
    axTxPSD = nexttile;
    grid(axTxPSD, 'on');
    hold(axTxPSD, 'on');
    xlabel(axTxPSD, 'Frequency (MHz)');
    ylabel(axTxPSD, 'Power (dBm)');
    title(axTxPSD, 'Transmitted Signal PSD');
    txPSDLine = plot(axTxPSD, NaN, NaN, 'b', 'LineWidth', 1.2);
    xlim(axTxPSD, [-Fs/2 Fs/2] / 1e6);

    %---------------- RECEIVED PSD AXIS ----------------
    axRxPSD = nexttile;
    grid(axRxPSD, 'on');
    hold(axRxPSD, 'on');
    xlabel(axRxPSD, 'Frequency (MHz)');
    ylabel(axRxPSD, 'Power (dBm)');
    title(axRxPSD, 'Received Signal PSD');
    rxPSDLine = plot(axRxPSD, NaN, NaN, 'r', 'LineWidth', 1.2);
    xlim(axRxPSD, [-Fs/2 Fs/2] / 1e6);

    % PSD settings based on the reference Comm_Sys_Baseband code.
    % The window is selected per block because this simulation uses
    % short blocks (12 symbols), shorter than 2048 samples.
    psdNFFT = 4096;
    psdResistance = 1;   % 1-ohm reference, as in the reference code

    %=========================================================
    %                  SIMULATION LOOP
    %=========================================================

    fprintf('Starting %s Transmission...\n', modScheme);
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

        % Random ACI data for every block
        aci_bits = randi([0 1], blockSize, 1);

        aci_syms = encoder(aci_bits, modScheme, M);

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
        % ACI block schedule
        %
        % [ACI_ON_blocks, ACI_OFF_blocks]
        %
        % During ACI-ON blocks, the complete random ACI
        % waveform is added to the desired signal.
        %
        % During ACI-OFF blocks, no ACI is added.
        %----------------------------------------------------

        blockPosition = mod(b - 1, aci_cycle_len);

        aciBlockOn = (blockPosition < aci_on_blocks);

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

        fprintf('Block %d: ACI block = %s\n', ...
            b, aciStatusText);

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

        noise = sqrt(noiseVar/2) * ...
            (randn(size(txAmp)) + 1j*randn(size(txAmp)));

        channelGain = 0.5;

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

        modScheme = adaptive_controller(isAck, 'immediate', false);

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
                refSyms = qammod(0:15, 16, ...
                    'UnitAveragePower', true);

            case '64-QAM'
                refSyms = qammod(0:63, 64, ...
                    'UnitAveragePower', true);

        end

        set(hIdeal, 'XData', real(refSyms), ...
            'YData', imag(refSyms));

        set(hScatter, 'XData', real(rxSyms), ...
            'YData', imag(rxSyms));

        maxVal = 1.5 * max(abs([real(refSyms(:)); ...
            imag(refSyms(:))]));

        axis(ax, [-maxVal maxVal -maxVal maxVal]);
        axis(ax, 'square');

        title(ax, sprintf('%s | Block %d/%d | Errors = %d | %s', ...
            currentModScheme, b, numBlocks, numErrors, status));

        %----------------------------------------------------
        % 10. UPDATE TRANSMITTED AND RECEIVED PSD
        %
        % Transmitted PSD: desired signal after TX amplifier,
        % BEFORE adding ACI.
        %
        % Received PSD: after channel, AWGN, RX amplifier and
        % RX nonlinearity; includes the effects of ACI.
        %----------------------------------------------------

        % Hann window, 50% overlap, 4096-point FFT, power spectrum.
        % Use the available block length because each block is short.
        psdWindowLength = min(2048, length(txAmp_noACI));
        psdWindow = hann(psdWindowLength, 'periodic');
        psdOverlap = floor(psdWindowLength / 2);

        % Transmitter PSD: desired signal only, before ACI addition
        [Ptx, fPSD] = pwelch(txAmp_noACI, psdWindow, psdOverlap, ...
            psdNFFT, Fs, 'centered', 'power');

        % Receiver PSD: received signal after channel and RX nonlinearity
        [Prx, ~] = pwelch(rxSignalForPSD, psdWindow, psdOverlap, ...
            psdNFFT, Fs, 'centered', 'power');

        % Convert power spectrum to dBm using a 1-ohm reference
        Ptx_dBm = 10*log10((Ptx / psdResistance) + eps) + 30;
        Prx_dBm = 10*log10((Prx / psdResistance) + eps) + 30;

        % Update both PSD plots
        set(txPSDLine, 'XData', fPSD/1e6, 'YData', Ptx_dBm);
        set(rxPSDLine, 'XData', fPSD/1e6, 'YData', Prx_dBm);

        title(axTxPSD, sprintf('Desired TX PSD (Before ACI) | Block %d', b));

        title(axRxPSD, sprintf('Received PSD (With ACI Effects) | Block %d', b));

        drawnow;
        pause(0.2);

        if ~isvalid(hFig)

            fprintf('Simulation stopped by user.\n');
            break;

        end

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

    ackRate = mean(ackHistory(1:completedBlocks));

    fprintf('\n========== FINAL SIMULATION RESULTS ==========\n');

    fprintf('Completed blocks        : %d\n', completedBlocks);

    fprintf('ACK rate                : %.2f%%\n', 100 * ackRate);

    fprintf('Total successful payload: %.0f bits\n', ...
        totalSuccessfulBits);

    fprintf('Average goodput         : %.4f bits/symbol\n', ...
        averageGoodput);

    fprintf('==============================================\n\n');

    %% =====================================================
    % MODULATION SCHEME VS BLOCK NUMBER
    % ======================================================

    schemeNames = {'BPSK', 'QPSK', '16-QAM', '64-QAM'};

    schemeLevels = zeros(completedBlocks, 1);

    for b = 1:completedBlocks

        schemeLevels(b) = find( ...
            strcmp(schemeNames, usedModScheme{b}), 1);

    end

    figure('Name', 'Adaptive Modulation History', ...
        'NumberTitle', 'off', 'Color', 'w');

    stairs(1:completedBlocks, schemeLevels, ...
        'b-', 'LineWidth', 1.8);

    hold on;

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

    title('Adaptive Modulation with ACK/NACK Feedback and ACI Schedule');

    yyaxis right;

    stairs(1:completedBlocks, ...
        double(aciHistory(1:completedBlocks)), ...
        'k--', 'LineWidth', 1.3);

    ylim([-0.1 1.1]);

    yticks([0 1]);
    yticklabels({'ACI OFF', 'ACI ON'});

    ylabel('ACI block state');

    yyaxis left;

    legend('Modulation used', 'NACK', 'ACI block state', ...
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

    persistent idx ackCount nackCount
    persistent histAck
    persistent blockCount lastSwitch
    persistent aimdScore

    persistent modHistory modObs probeCycle historyLocked
    persistent holdRemaining
    persistent scoreValue

    persistent gpHistory gpBlocksSinceProbe
    persistent gpProbeActive gpProbeBase gpProbeCount gpProbeAcks

    schemes = {'BPSK', 'QPSK', '16-QAM', '64-QAM'};

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

        modHistory = cell(1, 4);
        modObs = zeros(1, 4);

        for k = 1:4
            modHistory{k} = zeros(1, 0);
        end

        probeCycle = 0;
        historyLocked = false;
        holdRemaining = 0;

        scoreValue = 5;

        gpHistory = cell(1, 4);

        for k = 1:4
            gpHistory{k} = [];
        end

        gpBlocksSinceProbe = 0;

        gpProbeActive = false;
        gpProbeBase = 1;
        gpProbeCount = 0;
        gpProbeAcks = 0;

    end

    % Controller parameters

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

        %% -----------------------------------------------
        % IMMEDIATE ADAPTATION
        % -----------------------------------------------

        case 'immediate'

            if ack == 1
                idx = min(idx + 1, 4);
            else
                idx = max(idx - 1, 1);
            end

        %% -----------------------------------------------
        % THRESHOLD ADAPTATION
        % -----------------------------------------------

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

        %% -----------------------------------------------
        % PER-MODULATION ADAPTATION
        % -----------------------------------------------

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

        %% -----------------------------------------------
        % HISTORY-BASED ADAPTATION
        % -----------------------------------------------

        case 'history'

            modHistory{usedIdx} = [modHistory{usedIdx}, ack];

            if numel(modHistory{usedIdx}) > historyWindow

                modHistory{usedIdx} = ...
                    modHistory{usedIdx}(end-historyWindow+1:end);

            end

            modObs(usedIdx) = min( ...
                modObs(usedIdx) + 1, historyWindow);

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

        %% -----------------------------------------------
        % WINDOW-BASED ADAPTATION
        % -----------------------------------------------

        case 'window'

            histAck = [histAck ack];

            if numel(histAck) > Wwindow

                histAck = histAck(end-Wwindow+1:end);

            end

            if numel(histAck) >= Wwindow

                pACK = mean(histAck);

                canSwitch = ...
                    (blockCount - lastSwitch) >= dwellTime;

                newIdx = idx;

                if canSwitch && pACK >= TH

                    newIdx = min(idx + 1, 4);

                elseif canSwitch && pACK <= TL

                    newIdx = max(idx - 1, 1);

                end

                if newIdx ~= idx

                    idx = newIdx;

                    lastSwitch = blockCount;

                    histAck = [];

                end

            end

        %% -----------------------------------------------
        % AIMD ADAPTATION
        % -----------------------------------------------

        case 'aimd'

            if ack == 1

                aimdScore = min(aimdScore + alpha, xmax);

            else

                aimdScore = max(beta * aimdScore, xmin);

            end

            [~, idx] = min(abs(levels - aimdScore));

        %% -----------------------------------------------
        % SCORE-BASED ADAPTATION
        % -----------------------------------------------

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

        %% -----------------------------------------------
        % GOODPUT PROBE ADAPTATION
        % -----------------------------------------------

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
                            schemes{usedIdx}, probeGoodput, ...
                            baseGoodput, 100 * probeFailureRate);

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

                if gpBlocksSinceProbe >= probeInterval && ...
                        usedIdx < 4 && ~isempty(gpHistory{usedIdx})

                    gpProbeBase = usedIdx;

                    gpProbeActive = true;

                    gpProbeCount = 0;

                    gpProbeAcks = 0;

                    idx = usedIdx + 1;

                end

            end

        otherwise

            error(['Unknown method. Use immediate, threshold, ' ...
                'per_modulation, history, window, aimd, score, ' ...
                'or goodput_probe.']);

    end

    modScheme = schemes{idx};

end