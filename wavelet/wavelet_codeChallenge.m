clear

% signal is random noise
srate = 2048;
signal = smooth( randn(srate*6,1),3 );

clf
plot(signal)
plot(abs(fft(signal)))

%% FIR filter

% define filter parameters
lower_bnd = 10; % Hz
upper_bnd = 15; % Hz

tw = .1;

samprate  = 2048; % Hz
filtorder = 4*round(samprate/lower_bnd);

filter_shape = [ 0 0 1 1 0 0 ];
filter_freqs = [ 0 lower_bnd*(1-tw) lower_bnd ...
                 upper_bnd upper_bnd+upper_bnd*tw ...
                 (samprate/2) ] / (samprate/2);

filterkern = firls(filtorder,filter_freqs,filter_shape);

signalFIR = filtfilt(filterkern,1,signal);

%% wavelet

timevec = -1:1/srate:1;
freq = (lower_bnd+upper_bnd)/2;

csw = cos(2*pi*freq*timevec); % cosine wave
fwhm = .25; % full-width at half-maximum in seconds
gaussian = exp( -(4*log(2)*timevec.^2) / fwhm^2 ); % Gaussian

% Morlet wavelet
mw = csw .* gaussian/(2*pi*30);

signalMW = conv(signal,mw,'same');


% nc = length(mw)+length(signal)-1;
% kh = floor(length(mw)/2)+1;
% signalMW2 = ifft( fft(signal',nc).*fft(miwav,nc) );
% signalMW2 = signalMW2(kh:end-kh+1);



% plotting

%%

N = length(signal);

tv = (0:N-1)/srate;

clf
subplot(311)
plot(tv,signal)
xlabel('Time (s)')
title('Original signal')

subplot(312), cla; hold on
plot(tv,signalFIR)
plot(tv,signalMW)
xlabel('Time (s)')
title('Filtered signals')
legend({'FIR';'MW'})

subplot(313), cla, hold on
hz = linspace(0,srate,N);
plot(hz,abs(fft(signalFIR/N)))
plot(hz,abs(fft(signalMW/N)))
set(gca,'xlim',[0 20])
legend({'FIR';'MW'})
title('Amplitude spectra')
xlabel('Frequency (Hz)')

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% My Solution %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Goal: Use original signal and figure out to create
% an FIR filter and a morlet wavelength to generate filtered
% signals into amplitude spectra
% Apply knowledge with wavelet and filtering
clear
load wavelet_codeChallenge.mat
whos
pkg load signal

% Inside wavelet_codeChallenge.mat includes:
% signal: original signal
% signalFIR: filtered signal
% signalMW: signal used with morlet wavelet
% srate: sampling rate

% Display the original signal
figure;
subplot(4,1,1);
plot((0:length(signal)-1)/srate, signal);
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Parameters for FIR filter
lower_bnd = 10; % Lower bound frequency in Hz
upper_bnd = 15; % Upper bound frequency in Hz
tw = 0.1; % Transition width
samprate = srate; % Sampling rate
filtorder = 4 * round(samprate / lower_bnd); % Filter order

% Define filter shape and frequencies
filter_shape = [0 0 1 1 0 0];
filter_freqs = [0 lower_bnd * (1 - tw) lower_bnd ...
                upper_bnd upper_bnd + upper_bnd * tw ...
                (samprate / 2)] / (samprate / 2);

% Design FIR filter using least squares method
filterkern = firls(filtorder, filter_freqs, filter_shape);

% Apply FIR filter to the original signal
signalFIR = filtfilt(filterkern, 1, signal);

% Display the FIR filtered signal
subplot(4,1,2);
plot((0:length(signalFIR)-1)/srate, signalFIR);
title('FIR Filtered Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Parameters for Morlet wavelet
timevec = -1:1/srate:1; % Time vector
freq = (lower_bnd + upper_bnd) / 2; % Central frequency
fwhm = 0.25; % Full-width at half-maximum in seconds

% Generate Morlet wavelet
csw = cos(2 * pi * freq * timevec); % Cosine wave
gaussian = exp(-(4 * log(2) * timevec.^2) / fwhm^2); % Gaussian
mw = csw .* gaussian / (2 * pi * 30); % Morlet wavelet

% Apply Morlet wavelet to the original signal
signalMW = conv(signal, mw, 'same');

% Display the Morlet wavelet filtered signal
subplot(4,1,3);
plot((0:length(signalMW)-1)/srate, signalMW);
title('Morlet Wavelet Filtered Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Compute amplitude spectra
N = length(signal);
hz = linspace(0, srate/2, floor(N/2)+1);

% FFT of FIR filtered signal
fftFIR = fft(signalFIR) / N;
amp_spectrum_FIR = 2 * abs(fftFIR(1:floor(N/2)+1)); % Scale by 2 for single-sided spectrum

% FFT of Morlet wavelet filtered signal
fftMW = fft(signalMW) / N;
amp_spectrum_MW = 2 * abs(fftMW(1:floor(N/2)+1)); % Scale by 2 for single-sided spectrum

% Plot the amplitude spectra
subplot(4,1,4);
hold on;
plot(hz, amp_spectrum_FIR);
plot(hz, amp_spectrum_MW);
set(gca, 'xlim', [0 20]);
title('Amplitude Spectra');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
legend({'FIR', 'MW'});
hold off;
