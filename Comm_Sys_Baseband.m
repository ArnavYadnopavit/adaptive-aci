function Comm_Sys_Baseband()
    clc;
    close all;
    rng(10);

    % =========================================================
    % 1. SYSTEM PARAMETERS
    % =========================================================
    modScheme  = '64-QAM';
    numBlocks  = 20;
    numSymbols = 100;

    Rs        = 10e6;        % Symbol rate (10 MBaud)
    Bandwidth = 15e6;        % Target RF Bandwidth (15 MHz)
    sps       = 20;          % INCREASED sps to 16 (prevents aliasing/clipping)
    Fs        = Rs * sps;    % Baseband Sampling Frequency (160 MHz)

    % ACI parameters
    max_power       = 0.25;   % ACI/desired power ratio
    freq_separation = 15e6;  % Adjacent carrier separation = 15 MHz

    % Amplifier parameters
    ampGain_Tx        = 2;
    ampGain_Rx_linear = 1;

    % Baseband Equivalent Cubic nonlinear amplifier: y = a1*x + (3/4)*a3*|x|^2*x
    a1 = 1;
    a3 = 0.09;

    % Channel parameters
    channelGain = 0.5;
    noiseVar    = 1e-3;

    % =========================================================
    % 2. MODULATION SETUP & 3. RRC FILTER
    % =========================================================
    M = 64; % For 64-QAM
    k = log2(M);
    blockSize = numSymbols * k;

    span    = 6;
    rolloff = (Bandwidth / Rs) - 1;
    hRRC = rcosdesign(rolloff, span, sps, 'sqrt');
    hRRC = hRRC / norm(hRRC);

    % =========================================================
    % 4. CONSTELLATION FIGURE
    % =========================================================
    refSyms = qammod(0:M-1, M, 'UnitAveragePower', true);

    hFig = figure('Name', 'Live Constellation', 'NumberTitle', 'off', 'Color', 'w');
    ax = axes(hFig);
    hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'square');
    maxVal = max(abs([real(refSyms(:)); imag(refSyms(:))])) * 1.5;
    axis(ax, [-maxVal maxVal -maxVal maxVal]);
    xlabel(ax, 'In-Phase (I)'); ylabel(ax, 'Quadrature (Q)');
    title(ax, sprintf('Live %s Constellation', modScheme));

    plot(ax, real(refSyms), imag(refSyms), 'r*', 'MarkerSize', 10, 'DisplayName', 'Ideal');
    hScatter = plot(ax, NaN, NaN, 'b.', 'MarkerSize', 8, 'DisplayName', 'Received');
    legend(ax, 'Location', 'best');

    % Initialize accumulators for PSD averaging
    txAmp_noACI_all = [];
    txAmp_withACI_all = [];
    rxNL_noACI_all = [];
    rxNL_withACI_all = [];
    
    current_time = 0; % Track continuous time for ACI carrier

    % =========================================================
    % 5. SIMULATION LOOP
    % =========================================================
    fprintf('Starting %s transmission...\\n', modScheme);
    fprintf('Block\\tErrors\\tStatus\\n');
    fprintf('------------------------------\\n');

    for b = 1:numBlocks
        % -----------------------------------------------------
        % TRANSMITTER: DESIRED SIGNAL (Baseband)
        % -----------------------------------------------------
        txBits = randi([0 1], blockSize, 1);
        txSyms = encoder(txBits, modScheme, M);

        txWaveform = upfirdn(txSyms, hRRC, sps);
        txWaveform = txWaveform(:); % Complex baseband signal

        txBB = txWaveform;

        % -----------------------------------------------------
        % GENERATE ACI
        % -----------------------------------------------------
        aciBits = randi([0 1], blockSize, 1);
        aciSyms = encoder(aciBits, modScheme, M);

        aciWaveform = upfirdn(aciSyms, hRRC, sps);
        aciWaveform = aciWaveform(:);

        L = min(length(txBB), length(aciWaveform));
        txBB = txBB(1:L);
        aciWaveform = aciWaveform(1:L);
        
        % Create continuous time vector to avoid phase discontinuity splatter
        t = current_time + (0:L-1)' / Fs;
        current_time = current_time + L / Fs;

        % Shift ACI to +15 MHz natively in complex baseband
        aciBBraw = aciWaveform .* exp(1j * 2 * pi * freq_separation * t);

        desiredPower = mean(abs(txBB).^2);
        rawAciPower  = mean(abs(aciBBraw).^2);

        % Scale ACI
        aciBB = sqrt(max_power * desiredPower / max(rawAciPower, eps)) * aciBBraw;

        txBB_noACI   = txBB;
        txBB_withACI = txBB + aciBB;

        % -----------------------------------------------------
        % TRANSMITTER AMPLIFIER
        % -----------------------------------------------------
        txAmp_noACI   = ampGain_Tx * txBB_noACI;
        txAmp_withACI = ampGain_Tx * txBB_withACI;
        
        % Accumulate signals for smooth PSD calculation at the end
        txAmp_noACI_all   = [txAmp_noACI_all; txAmp_noACI];
        txAmp_withACI_all = [txAmp_withACI_all; txAmp_withACI];

        % -----------------------------------------------------
        % CHANNEL
        % -----------------------------------------------------
        % Complex Gaussian noise
        noise = sqrt(noiseVar/2) * (randn(size(txAmp_noACI)) + 1j*randn(size(txAmp_noACI)));

        rx_noACI   = channelGain * txAmp_noACI + noise;
        rx_withACI = channelGain * txAmp_withACI + noise;

        % -----------------------------------------------------
        % LINEAR RECEIVER AMPLIFIER
        % -----------------------------------------------------
        rxLinear_noACI   = ampGain_Rx_linear * rx_noACI;
        rxLinear_withACI = ampGain_Rx_linear * rx_withACI;

        % -----------------------------------------------------
        % NONLINEAR AMPLIFIER (Baseband Equivalent)
        % -----------------------------------------------------
        % Baseband cubic model: y = a1*x + (3/4)*a3*|x|^2*x
        rxNL_noACI   = a1 * rxLinear_noACI + (3/4) * a3 * (abs(rxLinear_noACI).^2) .* rxLinear_noACI;
        rxNL_withACI = a1 * rxLinear_withACI + (3/4) * a3 * (abs(rxLinear_withACI).^2) .* rxLinear_withACI;
        
        % Accumulate signals for smooth PSD calculation
        rxNL_noACI_all   = [rxNL_noACI_all; rxNL_noACI];
        rxNL_withACI_all = [rxNL_withACI_all; rxNL_withACI];

        % -----------------------------------------------------
        % DEMODULATION
        % -----------------------------------------------------
        rxBaseband = rxNL_withACI;
        rxFiltered = upfirdn(rxBaseband, hRRC, 1, sps);

        firstSymbol = span + 1;
        lastSymbol  = firstSymbol + numSymbols - 1;

        if length(rxFiltered) < lastSymbol
            continue;
        end

        rxSyms = rxFiltered(firstSymbol:lastSymbol);
        rxBits = decoder(rxSyms, modScheme, M);

        [isAck, numErrors] = generate_ack_nack(txBits, rxBits);
        if isAck
            status = 'ACK';
        else
            status = 'NACK';
        end

        fprintf('%d\\t\\t%d\\t\\t%s\\n', b, numErrors, status);

        set(hScatter, 'XData', real(rxSyms), 'YData', imag(rxSyms));
        title(ax, sprintf('%s | Block %d/%d | Errors = %d | %s', modScheme, b, numBlocks, numErrors, status));
        drawnow;
    end
    
    % Plot accumulated smooth PSDs at the end
    plot_signal_psd_comparison(txAmp_noACI_all, txAmp_withACI_all, Fs, '1. Transmitted Baseband PSD');
    plot_signal_psd_comparison(rxNL_noACI_all, rxNL_withACI_all, Fs, '2. Received PSD: Nonlinear Path');
end

function syms = encoder(bits, modScheme, M)
    k = log2(M);
    bitMatrix = reshape(bits, k, [])';
    intData = bi2de(bitMatrix, 'left-msb');
    syms = qammod(intData, M, 'UnitAveragePower', true);
end

function bits = decoder(syms, modScheme, M)
    rxInt = qamdemod(syms, M, 'UnitAveragePower', true);
    k = log2(M);
    bitMatrix = de2bi(rxInt, k, 'left-msb');
    bits = reshape(bitMatrix', [], 1);
end

function [isAck, errors] = generate_ack_nack(txBits, rxBits)
    errors = sum(txBits ~= rxBits);
    isAck = (errors == 0);
end

function plot_signal_psd_comparison(sigNoACI, sigWithACI, Fs, plotTitle)
    sigNoACI   = sigNoACI(:);
    sigWithACI = sigWithACI(:);

    windowLength = 2048; % Provides robust averaging for a smooth trace
    win          = hann(windowLength, 'periodic');
    noverlap     = floor(windowLength / 2);
    nfft         = 4096;

    [P1, f] = pwelch(sigNoACI, win, noverlap, nfft, Fs, 'centered');
    [P2, ~] = pwelch(sigWithACI, win, noverlap, nfft, Fs, 'centered');

    R = 50;
    P1_dBmHz = 10 * log10((P1 / R) + eps) + 30;
    P2_dBmHz = 10 * log10((P2 / R) + eps) + 30;

    figure('Color', 'w', 'Name', plotTitle);
    h1 = plot(f/1e6, P1_dBmHz, 'LineWidth', 1.5);
    hold on;
    h2 = plot(f/1e6, P2_dBmHz, 'LineWidth', 1.5);

    xlabel('Baseband Frequency (MHz)');
    ylabel('PSD (dBm/Hz)');
    title(plotTitle);
    legend([h1, h2], {'Without ACI', 'With ACI'}, 'Location', 'best');
    grid on;
    set(gca, 'FontSize', 11);
    
    % Note: axis is not clipped manually; it naturally spans [-Fs/2, Fs/2]
    
    hold off;
end