clc;
clear all;
close all;

% Define the binary data sequence
data = [1 0 1 1 0 0 1];

% Define the time vector
t = linspace(0, length(data), length(data) * 100);

% Unipolar NRZ
unipolar_nrz = repelem(data, 100);

% Polar NRZ
polar_nrz = repelem(2*data - 1, 100);

% Manchester
manchester = zeros(1, length(data) * 100);
for i = 1:length(data)
    if data(i) == 1
        manchester((i-1)*100+1:(i-1)*100+50) = 1;
        manchester((i-1)*100+51:i*100) = -1;
    else
        manchester((i-1)*100+1:(i-1)*100+50) = -1;
        manchester((i-1)*100+51:i*100) = 1;
    end
end

% Differential Manchester
diff_manchester = zeros(1, length(data) * 100);
last_level = -1;
for i = 1:length(data)
    if data(i) == 1
        diff_manchester((i-1)*100+1:(i-1)*100+50) = last_level;
        last_level = -last_level;
        diff_manchester((i-1)*100+51:i*100) = last_level;
    else
        diff_manchester((i-1)*100+1:(i-1)*100+50) = last_level;
        diff_manchester((i-1)*100+51:i*100) = -last_level;
    end
end

% AMI
ami = zeros(1, length(data)*100);
last_polarity = -1;
for i = 1:length(data)
    if data(i) == 1
        last_polarity = -last_polarity;
        ami((i-1)*100+1:i*100) = last_polarity;
    end
end

% Pseudoternary
pseudoternary = zeros(1, length(data)*100);
last_polarity = -1;
for i = 1:length(data)
    if data(i) == 0
        last_polarity = -last_polarity;
        pseudoternary((i-1)*100+1:i*100) = last_polarity;
    end
end

% Unipolar RZ
unipolar_rz = zeros(1, length(data) * 100);
for i = 1:length(data)
    if data(i) == 1
        unipolar_rz((i-1)*100+1:(i-1)*100+50) = 1;
    end
end

% Polar RZ
polar_rz = zeros(1, length(data) * 100);
for i = 1:length(data)
    if data(i) == 1
        polar_rz((i-1)*100+1:(i-1)*100+50) = 1;
    else
        polar_rz((i-1)*100+1:(i-1)*100+50) = -1;
    end
end

% Bipolar RZ
bipolar_rz = zeros(1, length(data) * 100);
last_polarity = -1;
for i = 1:length(data)
    if data(i) == 1
        last_polarity = -last_polarity;
        bipolar_rz((i-1)*100+1:(i-1)*100+50) = last_polarity;
    end
end

% Plot the waveforms
figure;
subplot(10,1,1);
plot(t, unipolar_nrz, 'LineWidth', 2);
title('Unipolar NRZ');
ylim([-1.5 1.5]);

subplot(10,1,2);
plot(t, polar_nrz, 'LineWidth', 2);
title('Polar NRZ');
ylim([-1.5 1.5]);

subplot(10,1,3);
plot(t, manchester, 'LineWidth', 2);
title('Manchester');
ylim([-1.5 1.5]);

subplot(10,1,4);
plot(t, ami, 'LineWidth', 2);
title('AMI');
ylim([-1.5 1.5]);

subplot(10,1,5);
plot(t, pseudoternary, 'LineWidth', 2);
title('Pseudoternary');
ylim([-1.5 1.5]);

subplot(10,1,6);
plot(t, diff_manchester, 'LineWidth', 2);
title('Differential Manchester');
ylim([-1.5 1.5]);

subplot(10,1,7);
plot(t, unipolar_rz, 'LineWidth', 2);
title('Unipolar RZ');
ylim([-1.5 1.5]);

subplot(10,1,8);
plot(t, polar_rz, 'LineWidth', 2);
title('Polar RZ');
ylim([-1.5 1.5]);

subplot(10,1,9);
plot(t, bipolar_rz, 'LineWidth', 2);
title('Bipolar RZ');
ylim([-1.5 1.5]);

% Function to find the maximum number of consecutive ones

% Function to find the maximum number of consecutive ones
function [max_ones,max_zeros] = find_max_consecutive_ones(signal)
    max_ones = 0;
    max_zeros = 0;
    current_ones = 0;
    current_zeros = 0;
    bit_interval = inf;
    last_transition = -1;

    % Find the bit interval
    count=0;
    bit_interval=100000;
    var=1;
    for i = 1:length(signal)-1
        if signal(i) == 1 && signal(i+1)==-1
            var=0;
        end
        if signal(i) == -1 && signal(i+1)==1
            var=1;
        end
        
        if var==0
            count=count+1;
        else
           
            if bit_interval>count && count~=0
                bit_interval=count;
                
            end
            count=0;
        end
    end
    
    
    
    % Reset variables for counting consecutive ones and zeros
    max_ones = 0;
    max_zeros = 0;
    current_ones = 0;
    current_zeros = 0;
    
    % Count consecutive ones and zeros using the bit interval
    for i = 1:bit_interval:length(signal)
        if signal(i) == 1
            current_ones = current_ones + 1;
            current_zeros = 0;
            if current_ones > max_ones
                max_ones = current_ones;
            end
        else
            current_zeros = current_zeros + 1;
            current_ones = 0;
            if current_zeros > max_zeros
                max_zeros = current_zeros;
            end
        end
    end
   
end

% Function to differentiate coding scheme
function coding_scheme = differentiate_coding(signal, epsilon, alpha, beta, papr)
    % Initialize variables
    level = 0;
    negative_values = false;
    zero_levels = false;

    % Determine the level for each sample
    for i = 1:length(signal)
        if signal(i) == 1 - epsilon
            level = 1;
        elseif signal(i) == 0 + epsilon
            level = 0;
        elseif signal(i) == -1
            level = -1;
        end

        % Update flags
        negative_values = negative_values || (level == -1);
        zero_levels = zero_levels || (level == 0);
    end
    
    % Determine coding scheme based on flags and PAPR
    [max_ones,max_zeros] = find_max_consecutive_ones(signal);
    if negative_values && zero_levels
        if papr == alpha
            coding_scheme = 'bipolar RZ';
        else
            if max_ones == 1
                if max_zeros==0
                    coding_scheme = 'polar RZ';
                else
                    coding_scheme = 'AMI or pseudoternary';
                end
            else
                coding_scheme = 'bipolar NRZ';
            end
        end
    elseif negative_values && ~zero_levels
        if max_ones <= 2
            coding_scheme = 'Manchester or Differential Manchester';
        else
            coding_scheme = 'polar NRZ';
        end
    elseif ~negative_values
        if papr > alpha
            coding_scheme = 'unipolar NRZ';
        else
            coding_scheme = 'unipolar RZ';
        end
    end
end


% Function to calculate PAPR
function papr = calculate_papr(signal)
    papr = max(signal.^2) / mean(signal.^2);
end

% Example usage
inputSignal = pseudoternaryZ; % Your input signal
PAPR = calculate_papr(inputSignal);
epsilon = 0; % Threshold for level decision
alpha = 0.25; % Threshold for PAPR comparison
beta = 2; % Threshold for consecutive 1's or 0's

coding_scheme = differentiate_coding(inputSignal, epsilon, alpha, beta, 1/PAPR);
disp(['Detected coding scheme: ', coding_scheme]);