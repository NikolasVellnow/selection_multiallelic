% Clear workspace
clear;
 
% Define variables
freqs = [0.005, 0.001, 0.01, 0.05, 0.1];
num_freqs = length(freqs);
num_x_values = 500;
 
% Pre-allocate a table for better performance
df = table('Size', [num_freqs * num_x_values, 3], ...
           'VariableTypes', {'double', 'double', 'double'}, ...
           'VariableNames', {'freq', 'coverage', 'prob_detection'});
 
counter = 0;
 
for i = 1:num_freqs
    freq = freqs(i);
    disp(freq);  % Equivalent to print in R
    for x = 1:num_x_values
        counter = counter + 1;
        df.freq(counter) = freq;
        df.coverage(counter) = x;
        prob = 1 - (1 - freq)^x;
        df.prob_detection(counter) = prob;
        disp(x);  % Equivalent to print in R
        disp(prob);
    end
end
 
% Convert 'freq' to a categorical variable
df.freq = categorical(df.freq);
 
% Plot the data
figure;
hold on;
colors = lines(num_freqs);  % Use different colors for each frequency
freqs = categorical(freqs);
for i = 1:num_freqs
    % Select data corresponding to the current frequency
    freq_data = df(df.freq == freqs(i), :);
    plot(freq_data.coverage, freq_data.prob_detection, 'Color', colors(i, :), 'Linewidth', 2);
end
hold off;
 
xlabel('Coverage');
ylabel('Probability of Detection');
legend(string(freqs), 'Location', 'best');
title('Probability of Detection vs Coverage');