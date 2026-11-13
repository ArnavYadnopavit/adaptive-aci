    
    function sweep_frequency_separation()
    
        %=========================================================
        % CONFIGURATION
        %=========================================================
    
        frequencySeparations = (0:0.25:16) * 1e6;
        numCases = numel(frequencySeparations);
    
        methods = {'history', 'q_learning', 'immediate'};
        numMethods = numel(methods);
    
        goodput = zeros(numCases, numMethods);
        ackRate = zeros(numCases, numMethods);
        successfulBits = zeros(numCases, numMethods);
    
        randomSeed = 12345;
    
        % Save original graphics settings
        oldVisibility = get(groot, 'DefaultFigureVisible');
    
        % Hide simulation figures
        set(groot, 'DefaultFigureVisible', 'off');
    
        % Ensure graphics settings are restored if an error occurs
        cleanupObj = onCleanup(@() set(groot, ...
            'DefaultFigureVisible', oldVisibility));
    
        fprintf('\n===== FREQUENCY-SEPARATION SWEEP =====\n');
    
        %=========================================================
        % FREQUENCY SEPARATION SWEEP
        %=========================================================
    
        for j = 1:numMethods
    
            fprintf('\nMethod: %s\n', methods{j});
    
            for i = 1:numCases
    
                fSep = frequencySeparations(i);
    
                % Same random seed for each method/frequency point
                rng(randomSeed, 'twister');
    
                % Store existing figures before simulation
                figuresBefore = findall(groot, 'Type', 'figure');
    
                % Run simulation
                evalc('[g,a,s] = Baseband_randomstr_with_aci_with_qlearning(fSep, methods{j});');
    
                % Store results
                goodput(i,j) = g;
                ackRate(i,j) = a;
                successfulBits(i,j) = s;
    
                fprintf('%.2f MHz | Goodput: %.4f | ACK: %.2f%%\n', ...
                    fSep/1e6, g, 100*a);
    
                % Close only figures created by this simulation
                figuresAfter = findall(groot, 'Type', 'figure');
                newFigures = setdiff(figuresAfter, figuresBefore);
    
                if ~isempty(newFigures)
                    close(newFigures);
                end
    
            end
        end
    
        %=========================================================
        % RESULTS TABLE
        %=========================================================
    
        results = table( ...
            frequencySeparations(:)/1e6, ...
            goodput(:,1), goodput(:,2), goodput(:,3), ...
            100*ackRate(:,1), 100*ackRate(:,2), 100*ackRate(:,3), ...
            'VariableNames', { ...
            'Frequency_MHz', ...
            'Goodput_History', ...
            'Goodput_QLearning', ...
            'Goodput_Immediate', ...
            'ACK_History', ...
            'ACK_QLearning', ...
            'ACK_Immediate'});
    
        fprintf('\n===== FINAL RESULTS =====\n');
        disp(results);
    
        %=========================================================
        % ENABLE GRAPHICS FOR FINAL PLOTS
        %=========================================================
    
        set(groot, 'DefaultFigureVisible', 'on');
    
        %=========================================================
        % GOODPUT COMPARISON
        %=========================================================
    
        figure('Name', 'Goodput vs Frequency Separation', ...
            'NumberTitle', 'off', 'Color', 'w');
    
        hold on;
        grid on;
    
        plot(frequencySeparations/1e6, goodput(:,1), ...
            '-o', 'Color', 'b', 'LineWidth', 1.8, ...
            'MarkerSize', 5, 'DisplayName', 'History');
    
        plot(frequencySeparations/1e6, goodput(:,2), ...
            '-s', 'Color', 'r', 'LineWidth', 1.8, ...
            'MarkerSize', 5, 'DisplayName', 'Q-Learning');
    
        plot(frequencySeparations/1e6, goodput(:,3), ...
            '-^', 'Color', 'g', 'LineWidth', 1.8, ...
            'MarkerSize', 5, 'DisplayName', 'Immediate');
    
        xlabel('Frequency Separation (MHz)');
        ylabel('Average Goodput (bits/symbol)');
        title('Goodput vs ACI Frequency Separation');
        legend('Location', 'best');
    
        hold off;
    
        %=========================================================
        % ACK RATE COMPARISON
        %=========================================================
    
        figure('Name', 'ACK Rate vs Frequency Separation', ...
            'NumberTitle', 'off', 'Color', 'w');
    
        hold on;
        grid on;
    
        plot(frequencySeparations/1e6, 100*ackRate(:,1), ...
            '-o', 'Color', 'b', 'LineWidth', 1.8, ...
            'MarkerSize', 5, 'DisplayName', 'History');
    
        plot(frequencySeparations/1e6, 100*ackRate(:,2), ...
            '-s', 'Color', 'r', 'LineWidth', 1.8, ...
            'MarkerSize', 5, 'DisplayName', 'Q-Learning');
    
        plot(frequencySeparations/1e6, 100*ackRate(:,3), ...
            '-^', 'Color', 'g', 'LineWidth', 1.8, ...
            'MarkerSize', 5, 'DisplayName', 'Immediate');
    
        xlabel('Frequency Separation (MHz)');
        ylabel('ACK Rate (%)');
        title('ACK Rate vs ACI Frequency Separation');
        legend('Location', 'best');
    
        hold off;
    
        %=========================================================
        % SAVE RESULTS
        %=========================================================
    
        writetable(results, 'frequency_separation_comparison.csv');
    
        save('frequency_separation_comparison.mat', ...
            'frequencySeparations', 'goodput', ...
            'ackRate', 'successfulBits', 'results', 'methods');
    
        fprintf('\nResults saved successfully.\n');
    
        % Restore original graphics setting
        clear cleanupObj;
    
    end
