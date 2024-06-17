%%
%     COURSE: Signal processing and image processing in MATLAB and Python
%    SECTION: Filtering
%      VIDEO: Code challenge: Filter these signals!
% Instructor: mikexcohen.com
%
%%

fs = 1000;
N = 10000;
x = randn(N,1);

%%% lowpass
fcutoff = 30;
transw  = .2;
order   = round( 5*fs/fcutoff );
shape   = [ 1 1 0 0 ];
frex    = [ 0 fcutoff fcutoff+fcutoff*transw fs/2 ] / (fs/2);


% filter
filtkern = firls(order,frex,shape);
y = filtfilt(filtkern,1,x);


%%% highpass
fcutoff = 5;
transw  = .05;
order   = round( 5*fs/fcutoff );
shape   = [ 0 0 1 1 ];
frex    = [ 0 fcutoff fcutoff+fcutoff*transw fs/2 ] / (fs/2);

% filter
filtkern = firls(order,frex,shape);
y = filtfilt(filtkern,1,y);


%%% notch
fcutoff = [ 18 24 ];
transw  = .1;
order   = round( 5*fs/fcutoff(1) );
shape   = [ 1 1 0 0 1 1 ];
frex    = [ 0 fcutoff(1)*(1-transw) fcutoff fcutoff(2)+fcutoff(2)*transw fs/2 ] / (fs/2);

% filter
filtkern = firls(order,frex,shape);
y = filtfilt(filtkern,1,y);

clf
subplot(211), hold on
plot(x,'r')
plot(y,'k')
title('Time domain')

yX = abs(fft(y)).^2;
xX = abs(fft(x)).^2;
hz = linspace(0,fs,N);

subplot(212), hold on
plot(hz,xX,'r');
plot(hz,yX,'k');
set(gca,'xlim',[0 80])
title('Frequency domain')

legend({'X';'Y'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% MY SOLUTION %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the data
load filtering_codeChallenge.mat
whos

% The data variables should be loaded now, fs, x, and y.
% fs: sampling frequency
% x: original signal (red in the plots)
% y: filtered signal (black in the plots)

%% Alternative Filter Design and Application using fir1

% Lowpass Filter using fir1
fcutoff = 30; % cutoff frequency in Hz
order   = round(5*fs/fcutoff); % filter order

% Design the lowpass filter
lowpass_filtkern = fir1(order, fcutoff/(fs/2), 'low');

% Apply the lowpass filter
y_low = filtfilt(lowpass_filtkern, 1, x);

% Highpass Filter using fir1
fcutoff = 5; % cutoff frequency in Hz
order   = round(5*fs/fcutoff); % filter order

% Design the highpass filter
highpass_filtkern = fir1(order, fcutoff/(fs/2), 'high');

% Apply the highpass filter
y_high = filtfilt(highpass_filtkern, 1, y_low);

% Notch Filter using fir1
fcutoff = [18 24]; % cutoff frequencies in Hz
order   = round(5*fs/fcutoff(1)); % filter order

% Design the notch filter
notch_filtkern = fir1(order, fcutoff/(fs/2), 'stop');

% Apply the notch filter
y_filtered = filtfilt(notch_filtkern, 1, y_high);

%% Plot the results in a new figure with multiple subplots

figure; % Create a new figure window

% Plot Time Domain
subplot(2, 1, 1), hold on;
plot(x, 'r');
plot(y_filtered, 'k');
title('Time Domain');
xlabel('Samples');
ylabel('Amplitude');
legend({'Original (X)', 'Filtered (Y)'});
grid on;

% Compute the Fourier Transforms
yX = abs(fft(y_filtered)).^2;
xX = abs(fft(x)).^2;
hz = linspace(0, fs, N);

% Plot Frequency Domain
subplot(2, 1, 2), hold on;
plot(hz, xX, 'r');
plot(hz, yX, 'k');
set(gca, 'xlim', [0 80]);
title('Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Power');
legend({'Original (X)', 'Filtered (Y)'});
grid on;

%%
