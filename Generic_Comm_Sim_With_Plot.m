function Generic_Comm_Sim_With_Plot()
    % --- System Parameters ---
    modScheme  = '16-QAM';   % Options: 'BPSK', 'QPSK', '16-QAM', '64-QAM'
    numBlocks  = 20;         % Increased to 20 to enjoy the animation
    numSymbols = 256;        % Number of SYMBOLS per block 
    snr_dB     = 12;         % Signal-to-Noise Ratio (dB)
    ampGain    = 1.5;        % Amplifier gain (linear)
    
    % --- Modulation Setup ---
    switch modScheme
        case 'BPSK',   M = 2;
        case 'QPSK',   M = 4;
        case '16-QAM', M = 16;
        case '64-QAM', M = 64;
        otherwise, error('Unsupported modulation scheme');
    end
    k = log2(M);             
    blockSize = numSymbols * k; 
    
    % --- RRC Filter Parameters ---
    sps     = 4;             
    span    = 6;             
    rolloff = 0.25;          
    hRRC    = rcosdesign(rolloff, span, sps, 'sqrt'); 

    % =========================================================
    % --- NEW: Visualization Setup ---
    % 1. Generate ideal reference constellation points
    switch modScheme
        case 'BPSK',   refSyms = pskmod(0:1, 2);
        case 'QPSK',   refSyms = pskmod(0:3, 4, pi/4);
        otherwise,     refSyms = qammod(0:M-1, M, 'UnitAveragePower', true);
    end
    
    % 2. Initialize the Figure
    hFig = figure('Name', 'Live Constellation Diagram', 'NumberTitle', 'off', 'Color', 'w');
    ax = axes(hFig);
    hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'square');
    
    % Scale axes slightly beyond the max constellation point for padding
    maxVal = max(abs([real(refSyms(:)); imag(refSyms(:))])) * 1.5;
    axis(ax, [-maxVal maxVal -maxVal maxVal]);
    
    title(ax, sprintf('Live %s Constellation (SNR = %d dB)', modScheme, snr_dB));
    xlabel(ax, 'In-Phase (I)');
    ylabel(ax, 'Quadrature (Q)');
    
    % Plot the ideal points once (Red Stars)
    plot(ax, real(refSyms), imag(refSyms), 'r*', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Ideal');
    
    % Prepare an empty plot handle for the noisy received points (Blue Dots)
    hScatter = plot(ax, NaN, NaN, 'b.', 'MarkerSize', 8, 'DisplayName', 'Received');
    legend(ax, 'Location', 'best');
    % =========================================================

    % --- Simulation Loop ---
    fprintf('Starting %s Transmission...\n', modScheme);
    fprintf('Block\tErrors\tStatus\n');
    fprintf('------------------------------\n');
    
    for b = 1:numBlocks
        %% 1. Transmitter
        txBits = randi([0 1], blockSize, 1);
        txSyms = encoder(txBits, modScheme, M);
        txWaveform = upfirdn(txSyms, hRRC, sps);
        
        %% 2. Channel & Amplifier
        txAmp = txWaveform * ampGain;
        rxWaveform = awgn(txAmp, snr_dB, 'measured');
        
        %% 3. Receiver
        rxFiltered = upfirdn(rxWaveform, hRRC, 1, sps);
        rxSyms = rxFiltered(span + 1 : span + numSymbols) / ampGain;
        rxBits = decoder(rxSyms, modScheme, M);
        
        %% 4. Metrics & Visualization
        [isAck, numErrors] = generate_ack_nack(txBits, rxBits);
        
        if isAck, status = 'ACK'; else, status = 'NACK'; end
        fprintf('%d\t\t%d\t\t%s\n', b, numErrors, status);
        
        % =========================================================
        % --- NEW: Update the scatterplot dynamically ---
        set(hScatter, 'XData', real(rxSyms), 'YData', imag(rxSyms));
        drawnow;        % Force MATLAB to refresh the UI immediately
        pause(0.2);     % Brief pause so you can actually see the animation
        % =========================================================
        
        % Break if the user closes the window early
        if ~isvalid(hFig)
            fprintf('Simulation stopped by user.\n');
            break;
        end
    end
end

%% --- Modular System Functions ---
function syms = encoder(bits, modScheme, M)
    k = log2(M);
    bitMatrix = reshape(bits, k, [])'; 
    intData = bi2de(bitMatrix, 'left-msb'); 
    
    switch modScheme
        case 'BPSK',   syms = pskmod(intData, 2);
        case 'QPSK',   syms = pskmod(intData, 4, pi/4); 
        otherwise,     syms = qammod(intData, M, 'UnitAveragePower', true);
    end
end

function bits = decoder(syms, modScheme, M)
    switch modScheme
        case 'BPSK',   rxInt = pskdemod(syms, 2);
        case 'QPSK',   rxInt = pskdemod(syms, 4, pi/4);
        otherwise,     rxInt = qamdemod(syms, M, 'UnitAveragePower', true);
    end
    k = log2(M);
    bitMatrix = de2bi(rxInt, k, 'left-msb');
    bits = reshape(bitMatrix', [], 1);
end

function [isAck, errors] = generate_ack_nack(txBits, rxBits)
    errors = sum(txBits ~= rxBits);
    if errors == 0, isAck = true; else, isAck = false; end
end