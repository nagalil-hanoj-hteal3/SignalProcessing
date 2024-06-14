%%
%     COURSE: Signal processing problems, solved in MATLAB and Python
%    SECTION: Time-domain denoising
%      VIDEO: Median filter to remove spike noise
% Instructor: sincxpress.com
%
%%

% create signal
n = 2000;
signal = cumsum(randn(n,1));


% proportion of time points to replace with noise
propnoise = .05;

% find noise points
noisepnts = randperm(n);
noisepnts = noisepnts(1:round(n*propnoise));

% generate signal and replace points with noise
signal(noisepnts) = 50+rand(size(noisepnts))*100;


% use hist to pick threshold
%figure(1), clf
%histogram(signal,100) <- matlab
%hist(signal,100) <- octave
%zoom on

% alternative for using the octave
% Use hist to pick threshold
figure(1); clf;
[counts, binCenters] = hist(signal, 100); % Get histogram counts and bin centers
bar(binCenters, counts); % Plot histogram as bar graph
zoom on;

% visual-picked threshold (important)
threshold = 40;

% find data values above the threshold
% Use This
suprathresh = find( signal>threshold );

%all of the timepoints should be replace median time points
%resulting to losing dynamics of the signal
%suprathresh = 1:n;

% initialize filtered signal
filtsig = signal;

% loop through suprathreshold points and set to median of k
k = 20; % actual window is k*2+1
for ti=1:length(suprathresh)

    % find lower and upper bounds
    lowbnd = max(1,suprathresh(ti)-k);
    uppbnd = min(suprathresh(ti)+k,n);

    % compute median of surrounding points
    filtsig(suprathresh(ti)) = median(signal(lowbnd:uppbnd));
end

% plot
figure(2), clf
plot(1:n,signal, 1:n,filtsig, 'linew',2)
zoom on

%% done.

