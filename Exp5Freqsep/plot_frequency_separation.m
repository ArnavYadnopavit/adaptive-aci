
function plot_frequency_separation()

    clc;
    close all;

    % Read previously saved results
    filename = 'frequency_separation_comparison.csv';

    if ~isfile(filename)
        error('CSV file not found: %s', filename);
    end

    results = readtable(filename);

    % Extract data
    fSep = results.Frequency_MHz;

    goodputHistory   = results.Goodput_History;
    goodputQLearning = results.Goodput_QLearning;
    goodputImmediate = results.Goodput_Immediate;

    ackHistory   = results.ACK_History;
    ackQLearning = results.ACK_QLearning;
    ackImmediate = results.ACK_Immediate;

    %=========================================================
    % GOODPUT COMPARISON
    %=========================================================

    figure('Name', 'Goodput vs Frequency Separation', ...
        'NumberTitle', 'off', 'Color', 'w');

    plot(fSep, goodputHistory, '-o', ...
        'Color', 'b', 'LineWidth', 1.8, ...
        'MarkerSize', 5, 'DisplayName', 'History');

    hold on;

    plot(fSep, goodputQLearning, '-s', ...
        'Color', 'r', 'LineWidth', 1.8, ...
        'MarkerSize', 5, 'DisplayName', 'Q-Learning');

    plot(fSep, goodputImmediate, '-^', ...
        'Color', 'g', 'LineWidth', 1.8, ...
        'MarkerSize', 5, 'DisplayName', 'Immediate');

    grid on;
    xlabel('Frequency Separation (MHz)');
    ylabel('Average Goodput (bits/symbol)');
    title('Goodput vs ACI Frequency Separation');

    legend('Location', 'best');
    xlim([min(fSep), max(fSep)]);

    hold off;

    %=========================================================
    % ACK RATE COMPARISON
    %=========================================================

    figure('Name', 'ACK Rate vs Frequency Separation', ...
        'NumberTitle', 'off', 'Color', 'w');

    plot(fSep, ackHistory, '-o', ...
        'Color', 'b', 'LineWidth', 1.8, ...
        'MarkerSize', 5, 'DisplayName', 'History');

    hold on;

    plot(fSep, ackQLearning, '-s', ...
        'Color', 'r', 'LineWidth', 1.8, ...
        'MarkerSize', 5, 'DisplayName', 'Q-Learning');

    plot(fSep, ackImmediate, '-^', ...
        'Color', 'g', 'LineWidth', 1.8, ...
        'MarkerSize', 5, 'DisplayName', 'Immediate');

    grid on;
    xlabel('Frequency Separation (MHz)');
    ylabel('ACK Rate (%)');
    title('ACK Rate vs ACI Frequency Separation');

    legend('Location', 'best');
    xlim([min(fSep), max(fSep)]);

    hold off;

    %=========================================================
    % PRINT SUMMARY
    %=========================================================

    fprintf('\n===== CSV RESULTS SUMMARY =====\n');

    fprintf('Number of frequency points: %d\n', height(results));

    fprintf('\nMaximum Goodput:\n');
    fprintf('History    : %.4f\n', max(goodputHistory));
    fprintf('Q-Learning : %.4f\n', max(goodputQLearning));
    fprintf('Immediate  : %.4f\n', max(goodputImmediate));

    fprintf('\nMaximum ACK Rate:\n');
    fprintf('History    : %.2f%%\n', max(ackHistory));
    fprintf('Q-Learning : %.2f%%\n', max(ackQLearning));
    fprintf('Immediate  : %.2f%%\n', max(ackImmediate));

    fprintf('\nPlots generated from saved CSV successfully.\n');

end
