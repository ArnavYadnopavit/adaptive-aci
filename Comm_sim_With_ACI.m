function Generic_Comm_Sim_With_Plot()

    % --- System Parameters ---
    modScheme  = '64-QAM';   % 'BPSK', 'QPSK', '16-QAM', '64-QAM'
    numBlocks  = 20;
    numSymbols = 10;
    fc = 2.4e9;   % Carrier frequency
    Fs = 50e9;    % Sampling frequency
    Rs = 1e6;     % Symbol rate

    %ACI parameters
    max_power = 0.05 ;
    freq_seperation = 1.25e6/2;
    aci_bandwidth = 1.25e6;
    max_duty_cycle = 0.4 ;
    adjacentCenter = fc + freq_seperation;

    % Amplifiers
    ampGain_Tx          = 2;
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




    % =========================================================
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
    plot(ax, ...
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


    for b = 1:numBlocks

        %% ---------------------------------------------------
        % 1. TRANSMITTER
        % ---------------------------------------------------

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

        % Randomize the interference with a maximum duty cycle

        % Calculate maximum allowed interfered symbols
        maxInterfered = floor(max_duty_cycle * numSymbols);

        % Initially, no symbols are interfered
        activeSymbols = false(numSymbols, 1);

        % Randomly choose the number of interfered symbols (0 to maxInterfered)
        numInterfered = randi([0, maxInterfered]);

        % Randomly select distinct symbol positions
        if numInterfered > 0
            selected = randperm(numSymbols, numInterfered);
            activeSymbols(selected) = true;
        end

        % Create a sample-level mask
        N = numel(aci_RF);
        sampleIndex = (0:N-1).';

        % Map samples to symbol intervals
        symbolIndex = min(floor(sampleIndex / sps) + 1, numSymbols);

        % Apply the symbol selections to the waveform
        aciMask = double(activeSymbols(symbolIndex));

        % Apply interference mask
        aci_interference = aciMask .* aci_RF;

        % Add interference to the desired RF signal
        txRF = txRF + aci_interference;

        % Display the number of interfered symbols
        fprintf('Interfered symbols: %d out of %d\n', ...
            sum(activeSymbols), numSymbols);
        
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


        if isAck
            status = 'ACK';
        else
            status = 'NACK';
        end


        fprintf('%d\t\t%d\t\t%s\n', ...
            b, numErrors, status);


        %% ---------------------------------------------------
        % 9. UPDATE LIVE CONSTELLATION
        % ---------------------------------------------------

        set(hScatter, ...
            'XData', real(rxSyms), ...
            'YData', imag(rxSyms));


        title(ax, sprintf( ...
            '%s | Block %d/%d | Errors = %d | %s', ...
            modScheme, ...
            b, ...
            numBlocks, ...
            numErrors, ...
            status));


        drawnow;


        % Small delay so animation can be seen
        pause(0.2);


        % Stop if figure is closed
        if ~isvalid(hFig)
            fprintf('Simulation stopped by user.\n');
            break;
        end

    end

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
