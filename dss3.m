clear all; close all; clc;

%% Parameters
Fs = 1e6;               % Sampling frequency (1 MHz)
fc = 1e5;               % Carrier frequency (100 kHz)
bit_t = 1e-3;           % Bit duration (1 ms)
num_bits = 100;         % Number of bits in the pseudo-random code
SNR = 10;               % Signal-to-Noise Ratio in dB
num_targets = 3;        % Number of targets

%% Generate the Pseudo-Random Code (PN code) for DSSS
pn_code = randi([0, 1], 1, num_bits);  % Generate random bits
pn_code(pn_code == 0) = -1;            % Map 0s to -1 for BPSK modulation

%% Message Signal (DSSS) - Example message generation
message = [0 1 1 0 1 0 1 1];           % Example binary message
message(message == 0) = -1;            % Map 0s to -1 for BPSK modulation

% Adjust num_bits to ensure it's a multiple of message length
message_length = length(message);
if mod(num_bits, message_length) ~= 0
    num_bits = message_length * ceil(num_bits / message_length); % Adjust num_bits to the nearest multiple
end

% Expand the message to match the length required for spreading
message_signal = kron(message, ones(1, num_bits / message_length));

% Ensure pn_code matches the length of message_signal
if length(pn_code) < length(message_signal)
    pn_code = repmat(pn_code, 1, ceil(length(message_signal) / length(pn_code)));
end
pn_code = pn_code(1:length(message_signal));  % Truncate if necessary

spread_signal = message_signal .* pn_code;  % Spread the message using the PN code

%% Carrier Signal Generation
t = (0:length(spread_signal)-1) / Fs;  % Time vector
carrier = cos(2 * pi * fc * t);        % Carrier signal

%% Transmit the Spread Spectrum Signal
tx_signal = spread_signal .* carrier;  % Modulated transmit signal

%% Simulate the Received Signal with Multiple Target Reflections and Noise
delays = [50, 150, 300];  % Delays in samples representing the targets' distances
rx_signal = zeros(1, length(tx_signal));

% Simulate reflections from multiple targets
for i = 1:num_targets
    % Create the delayed signal by padding zeros at the beginning
    delayed_signal = [zeros(1, delays(i)), tx_signal];
    
    % Ensure the delayed signal has the same length as rx_signal
    delayed_signal = delayed_signal(1:length(rx_signal));
    
    % Add the delayed signal to the received signal
    rx_signal = rx_signal + delayed_signal;  % Accumulate reflections from each target
end

%% Add Noise to the Received Signal
% Calculate the power of the transmitted signal
signal_power = mean(abs(tx_signal).^2);

% Calculate the noise power based on the desired SNR
noise_power = signal_power / (10^(SNR/10));

% Generate Gaussian noise with the calculated noise power
noise = sqrt(noise_power) * randn(size(tx_signal));

% Add the noise to the received signal
noisy_signal = rx_signal + noise;

%% Despreading and Demodulation using a Matched Filter
matched_filter = fliplr(pn_code);  % Matched filter is the time-reversed PN code
demod_signal = noisy_signal .* carrier;  % Demodulate the received signal

% Perform matched filtering on the demodulated signal
matched_output = conv(demod_signal, matched_filter, 'same');

%% Original Code: DSSS Modulation (Retain Original Plots)
Fs_old = 1000; fc_old = 100; fp = 4; bit_t_old = 0.1;
m = [0 0 1 1 1 1 0 0];  % Original message
m(m == 0) = -1;         % BPSK mapping
message_old = repmat(m, fp, 1);
message_old = reshape(message_old, 1, []);

pn_code_old = randi([0, 1], 1, length(m) * fp);
pn_code_old(pn_code_old == 0) = -1;
DSSS_old = message_old .* pn_code_old;

%% Plot Time-Domain Signal, PN Code, and DSSS Signal from Original Code
figure;
subplot(3,1,1); stairs(message_old, 'linewidth', 2);
title('Message bit sequence (Original)');
axis([0 length(message_old) -1 1]);

subplot(3,1,2); stairs(pn_code_old, 'linewidth', 2);
title('Pseudo-random code (Original)');
axis([0 length(pn_code_old) -1 1]);

subplot(3,1,3); stairs(DSSS_old, 'linewidth', 2);
title('Modulated signal (Original)');
axis([0 length(DSSS_old) -1 1]);

%% Plot Frequency-Domain Signals from Original Code
f_old = linspace(-Fs_old/2, Fs_old/2, 1024);
figure;
subplot(3,1,1); plot(f_old, abs(fftshift(fft(message_old, 1024))), 'linewidth', 2);
title('Message spectrum (Original)');

subplot(3,1,2); plot(f_old, abs(fftshift(fft(pn_code_old, 1024))), 'linewidth', 2);
title('Pseudo-random code spectrum (Original)');

subplot(3,1,3); plot(f_old, abs(fftshift(fft(DSSS_old, 1024))), 'linewidth', 2);
title('Modulated signal spectrum (Original)');

%% DSSS Modulation Figures (New Code)
figure;
subplot(3, 1, 1);
stairs(message_signal, 'linewidth', 2);
title('Original Message Signal (New)');
xlabel('Sample Index'); ylabel('Amplitude');
axis([0 length(message_signal) -1.5 1.5]);

subplot(3, 1, 2);
stairs(pn_code, 'linewidth', 2);
title('Pseudo-Random Code (PN Code) (New)');
xlabel('Sample Index'); ylabel('Amplitude');
axis([0 length(pn_code) -1.5 1.5]);

subplot(3, 1, 3);
stairs(spread_signal, 'linewidth', 2);
title('DSSS Spread Signal (New)');
xlabel('Sample Index'); ylabel('Amplitude');
axis([0 length(spread_signal) -1.5 1.5]);

%% Matched Filter and Radar Target Detection Figures (New Code)
figure;
subplot(3, 1, 1);
plot(tx_signal);
title('Transmitted DSSS Signal (Spread Spectrum) (New)');
xlabel('Sample Index'); ylabel('Amplitude');

subplot(3, 1, 2);
plot(noisy_signal);
title('Received Signal with Multiple Target Reflections and Noise (New)');
xlabel('Sample Index'); ylabel('Amplitude');

subplot(3, 1, 3);
plot(abs(matched_output));
title('Matched Filter Output for Target Detection (New)');
xlabel('Sample Index'); ylabel('Amplitude');
hold on;
threshold = 0.6 * max(abs(matched_output));  % Threshold for target detection
plot([1 length(matched_output)], [threshold threshold], 'r--');
legend('Matched Filter Output', 'Detection Threshold');

%% Identify Targets Based on Threshold (New Code)
detected_targets = find(abs(matched_output) > threshold);
disp('Detected Targets at Sample Indices:');
disp(detected_targets);

% Visualize target positions
figure;
stem(detected_targets, ones(size(detected_targets)), 'r', 'LineWidth', 2);
title('Detected Targets at Different Ranges');
xlabel('Sample Index'); ylabel('Target Presence');
