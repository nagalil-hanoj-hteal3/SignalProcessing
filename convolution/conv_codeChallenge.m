% conv first with planck then wavelet in there, vs. wavelet only.
% different? why or why not, then test with data

%% Goal: Reproduce the mean filter using convolution during frequency domain
% of this method
% k = 20;
% for i = k+1 : n-k-1
%  filtsig(i) = mean(signal(i-k:i+k)); % <- each point is average of k surrounding points
% end
% Come up with convolution kernel and take the fourier spectrum
% of the kernel times the fourier spectrum of the signal and
% take the inverse of fourier transform of the spectral multiplication
% that should be almost the same as the method

% Generate broadband noise
N = 10000;
signal = randn(N, 1);

% Define the mean filter kernel
k = 20;
kernel = ones(2*k+1, 1) / (2*k+1);

% Zero-pad the kernel to the length of the signal for frequency domain convolution
kernel_padded = [kernel; zeros(N - length(kernel), 1)];

% Fourier transform of the signal and the kernel
signal_fft = fft(signal);
kernel_fft = fft(kernel_padded);

% Perform spectral multiplication
filtered_signal_fft = signal_fft .* kernel_fft;

% Inverse Fourier transform to get the filtered signal in time domain
filtered_signal = ifft(filtered_signal_fft);

% Compare with the time-domain mean filter
filtsig = signal; % Initialize the filtered signal
for i = k+1 : N-k
    filtsig(i) = mean(signal(i-k:i+k));
end

% Plot both signals to compare
figure;
subplot(2, 1, 1);
plot(real(filtered_signal));
title('Filtered Signal (Frequency Domain Convolution)');
subplot(2, 1, 2);
plot(filtsig);
title('Filtered Signal (Time Domain Mean Filter)');

% Plot both signals in the same figure
%figure;
%plot(real(filtered_signal), 'b');
%hold on;
%plot(filtsig, 'r');
%title('Comparison of Filtered Signals');
%legend('Frequency Domain Convolution', 'Time Domain Mean Filter');
%xlabel('Sample Index');
%ylabel('Amplitude');
%hold off;
