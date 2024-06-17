%%
%     COURSE: Signal processing problems, solved in MATLAB and Python
%    SECTION: Filtering
%      VIDEO: Use filtering to separate birds in a recording
% Instructor: sincxpress.com
%
%%

%% Clear and setup

%%%%% MOST OF THE CODE REVISED INTO OCTAVE GNU NOT MATLAB IDE

clc;
clear;

% Ensure the signal package is installed and loaded in Octave
pkg load signal

% Read the audio file
[bc, fs] = audioread('XC403881.wav');
N = length(bc);

%% Compute and plot the spectrogram
window_size = 1000;
noverlap = 100;
nfft = 1000;

window = hann(window_size);
num_windows = floor((N - noverlap) / (window_size - noverlap));
powspect = zeros(nfft/2+1, num_windows);

for i = 1:num_windows
    start_index = (i-1) * (window_size - noverlap) + 1;
    end_index = start_index + window_size - 1;
    segment = bc(start_index:end_index, 2) .* window;
    fft_segment = fft(segment, nfft);
    powspect(:, i) = abs(fft_segment(1:nfft/2+1)).^2;
end

time = (0:num_windows-1) * (window_size - noverlap) / fs;
frex = (0:nfft/2) * fs / nfft;

figure(1), clf;
imagesc(time, frex, powspect);
axis xy;
set(gca, 'clim', [0 1]*2, 'ylim', [0 15000], 'xlim', [0 max(time)]);
xlabel('Time (sec.)'), ylabel('Frequency (Hz)');
colormap hot;
colorbar;

%% Select frequency ranges based on visual inspection
frange{1} = [1700 2600];
frange{2} = [5100 6100];

% Draw boundary lines on the plot
colorz = 'wm';
hold on;
for fi = 1:length(frange)
    plot(get(gca,'xlim'), [1 1] * frange{fi}(1), [colorz(fi) '--']);
    plot(get(gca,'xlim'), [1 1] * frange{fi}(2), [colorz(fi) '--']);
end

%% Compute and apply FIR filters
filteredSig = cell(2,1);

for filteri = 1:length(frange)
    order = round(10 * fs / frange{1}(1));
    filtkern = fir1(order, frange{filteri} / (fs / 2));

    for chani = 1:2
        dat1chan = bc(:, chani);

        % Zero-phase-shift filter with reflection
        sigR = [dat1chan(end:-1:1); dat1chan; dat1chan(end:-1:1)]; % reflect
        fsig = filter(filtkern, 1, sigR); % forward filter
        fsig = filter(filtkern, 1, fsig(end:-1:1)); % reverse filter
        fsig = fsig(end:-1:1); % reverse again for 0 phase
        fsig = fsig(N+1:end-N); % chop off reflected parts

        filteredSig{filteri}(:, chani) = fsig;
    end
end

%% Play the original and filtered signals
disp('Playing original signal...');
soundsc(bc, fs);

%pause(length(bc) / fs + 2); % Wait for original sound to finish

disp('Playing lower frequency range...');
soundsc(filteredSig{1}, fs);

%pause(length(filteredSig{1}) / fs + 2); % Wait for first filtered sound to finish

disp('Playing higher frequency range...');
soundsc(filteredSig{2}, fs);
