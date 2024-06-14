%%
% COURSE: Signal processing and image processing in MATLAB and Python
% SECTION: Time-domain denoising
% VIDEO: Code challenge: Denoise these signals!
% Instructor: mikexcohen.com
%
%%

%% added
clc;
clear;

% Load the signal package
pkg load signal;

N = 4000;

origSignal = linspace(-1,1,N) .* sin(linspace(0,10*pi,N)) + randn(1,N);
r = randperm(N);
nn = round(N*.05);
origSignal(r(1:nn)) = (1+rand(1,nn))*10;
origSignal(r(end-nn+1:end)) = -(1+rand(1,nn))*10;

% Initialize cleanedSignal with origSignal
cleanedSignal = origSignal;

% remove positive noise spikes
p2r = find(origSignal > 5);
k = 5;
for i = 1:length(p2r)
    cleanedSignal(p2r(i)) = median(cleanedSignal(max(1, p2r(i) - k):min(N, p2r(i) + k)));
end

% remove negative noise spikes
p2r = find(origSignal < -5);
k = 5;
for i = 1:length(p2r)
    cleanedSignal(p2r(i)) = median(cleanedSignal(max(1, p2r(i) - k):min(N, p2r(i) + k)));
end

% mean-smooth
k = 150;
for i = 1:N
    cleanedSignal(i) = mean(cleanedSignal(max(1, i - k):min(N, i + k)));
end

% Ensure the figure is cleared before plotting
figure;
subplot(3,1,1);
plot(1:N, origSignal, 'linew', 2);
title("Original Signal");

subplot(3,1,2);
plot(1:N, cleanedSignal, 'linew', 2);
title("Lecturer's Solution");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Your Result: %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load your saved data
load denoising_codeChallenge.mat;

% Ensure the size matches
N = length(origSignal);

% Step 1: Identify and replace outliers
k = 5; % Window size for median filter
outlierThreshold = 5; % Threshold for identifying outliers
myCleanedSignal = origSignal;

% Positive outliers
outlierIndices = find(origSignal > outlierThreshold);
for i = 1:length(outlierIndices)
    windowStart = max(1, outlierIndices(i) - k);
    windowEnd = min(N, outlierIndices(i) + k);
    myCleanedSignal(outlierIndices(i)) = median(origSignal(windowStart:windowEnd));
end

% Negative outliers
outlierIndices = find(origSignal < -outlierThreshold);
for i = 1:length(outlierIndices)
    windowStart = max(1, outlierIndices(i) - k);
    windowEnd = min(N, outlierIndices(i) + k);
    myCleanedSignal(outlierIndices(i)) = median(origSignal(windowStart:windowEnd));
end

% Step 2: Moving average smoothing
windowSize = 150; % Window size for moving average (odd number)
movingAvgFilter = ones(1, windowSize) / windowSize; % Normalized moving average filter
mySmoothSignal = conv(myCleanedSignal, movingAvgFilter, 'same'); % Apply moving average filter

subplot(3,1,3);
plot(1:N, mySmoothSignal, 'Color', [0.8 0.2 0.1], 'linew', 2);
title('Your Solution');

%% done.

