clc;

% create signal
srate = 1000;
time  = -3:1/srate:3;
pnts  = length(time);
freqmod = exp(-time.^2)*10+10;
freqmod = freqmod + linspace(0,10,pnts);
signal  = sin( 2*pi * ((time + cumsum(freqmod))/srate) );


% plot the signal
figure(1), clf
subplot(3,1,1)
plot(time,signal,'linew',1)
xlabel('Time (s)')
title('Time-domain signal')


n  = 500;
hz = linspace(0,srate,n+1);
tf = zeros(floor(pnts/n)-1,length(hz));
tv = zeros(floor(pnts/n)-1,1);

for i=1:floor(pnts/n)-1
    datasnip = signal(i*n:(i+1)*n);

    pw = abs(fft(datasnip)).^2;
    tf(i,1:length(hz)) = pw(1:length(hz));
    tv(i) = mean(time(i*n:(i+1)*n));
end

subplot(3,1,2)
imagesc(tv,hz,tf')
axis xy
set(gca,'ylim',[0 40])
colormap hot
xlabel('Time (s)'), ylabel('Frequency (Hz)')
title("Professor's Solution")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% MY SOLUTION %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load spectral_codeChallenge.mat;

n = length(signal);
winlength = srate / 2;
td_wlen_start = 1:winlength:n-winlength; % Time start of each window

hzW = linspace(0, srate/2, floor(winlength));
signalpowW = [];

% hannw = 0.5 - cos(2*pi*linspace(0,1,winlength))./2;
for wi = 1:length(td_wlen_start)
    datachunk = signal(td_wlen_start(wi):td_wlen_start(wi)+winlength-1);
    % datachunk = datachunk.*hannw;
    tmppow = abs(fft(datachunk)/winlength).^2;
    signalpowW = [signalpowW, tmppow'];
end

% Create a time vector based on the window start times
time_vector = (td_wlen_start - 1) / srate;

subplot(3,1,3);
imagesc(time_vector, hzW, signalpowW);
axis xy;
ylim([0 40]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('My Solution');
colorbar;
