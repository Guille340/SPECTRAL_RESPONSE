% Input
sourceName = 'AKP_3400LLX_080_2000_80';
sourceNameTitle = 'AKP 3400LLX 080 2000 80';
fileName = 'AKP_3400LLX_080_2000_80 (fs = 2kHz).csv';

% Extract Signal
csvData = csvread(fileName,1,0);
t = csvData(:,1)*1e-3; % time [s]
x_gun = csvData(:,2); % pressure [bar m]
nSamples = length(x_gun);
sampleInterval = unique(round(diff(t)*1e6)*1e-6);
win = hanning(nSamples);

% General Parameters
pref = 1e-6;
bar2pascal = 1e5;
fs = 1/sampleInterval;

% Sinewave
f = [100.09765625 250];
x_sin = sin(2*pi*f(1)*t) + 0.5*sin(2*pi*f(2)*t);

% Wrapping Effect
X = fft(x_gun); % original FFT
f = (0:length(X)-1)*fs/length(X);
x_wra256 = datawrap(x_gun,256);
X_wra256 = fft(x_wra256); % wrapped FFT
f_wra256 = (0:length(X_wra256)-1)*fs/length(X_wra256);
X_tru256 = fft(x_gun,256); % truncated FFT
f_tru256 = (0:length(X_tru256)-1)*fs/length(X_tru256);

figure
hold on
plot(f,20*log10(abs(X)),'LineWidth',1.5,'Color','b')
plot(f_wra256,20*log10(abs(X_wra256)),'LineWidth',1.5,'Color','m')
plot(f_tru256,20*log10(abs(X_tru256)),'LineWidth',1.5,'Color','g')
xlim([0 250])
xlabel('Frequency [Hz]')
ylabel('Relative Level [dB]')
legend({'Original (N = L = 4001)','Datawrap (N = 256)',...
    'Truncated (N = 256)'},'Location','SouthEast')
box on
set(gcf,'PaperPositionMode','auto')
title('Wrapping Effect')
print('Wrapping Effect','-dpng','-r250')
savefig('Wrapping Effect')

% Zero-Padding Effect
X = fft(x_gun); % original FFT
f = (0:length(X)-1)*fs/length(X);
X_pad = fft(x_gun,16384); % padded FFT
f_pad = (0:length(X_pad)-1)*fs/length(X_pad);

figure
hold on
plot(f,20*log10(abs(X)),'LineWidth',1.5,'Color','b')
plot(f_pad,20*log10(abs(X_pad)),'LineWidth',1,'Color','m')
xlim([0 250])
xlabel('Frequency [Hz]')
ylabel('Relative Level [dB]')
legend({'Original (N = L = 4001)','Zero-Pad (N = 16384)'},...
    'Location','SouthEast')
box on
set(gcf,'PaperPositionMode','auto')
title('Zero-Padding Effect')
print('Zero-Padding Effect','-dpng','-r250')
savefig('Zero-Padding Effect')

figure
hold on
plot(f,20*log10(abs(X)),'LineWidth',1.5,'Color','b')
plot(f_pad,20*log10(abs(X_pad)),'LineWidth',1,'Color','m')
xlim([5 20])
ylim([45 60])
xlabel('Frequency [Hz]')
ylabel('Relative Level [dB]')
legend({'Original (N = L = 4001)','Zero-Pad (N = 16384)'},...
    'Location','SouthEast')
box on
set(gcf,'PaperPositionMode','auto')
title('Zero-Padding Effect (Detail)')
print('Zero-Padding Effect (Detail)','-dpng','-r250')
savefig('Zero-Padding Effect (Detail)')

% Window Gain Effect (Sine)
nfft = 8192;
x_win = x_sin.*win;
winGain1 = rms(win);
winGain2 = mean(win);
X = fft(x_win,nfft)/winGain2; % original FFT
X = X.*conj(X)/nSamples^2;
X = [X(1); 2*X(2:ceil(nfft/2))];
f = (0:length(X)-1)*fs/(2*length(X));

figure
hold on
plot(f,10*log10(X*winGain2/winGain1),'LineWidth',1,'Color','b')
plot(f,10*log10(X),'LineWidth',1,'Color','m')
xlabel('Frequency [Hz]')
ylabel('Power Spectrum [dB re 1V]')
legend({'W_{NG}','W_{CG}'},...
    'Location','NorthEast')
box on
set(gcf,'PaperPositionMode','auto')
title('Window Gain Effect (CG)')
print('Window Gain Effect (CG)','-dpng','-r250')
savefig('Window Gain Effect (CG)')

figure
hold on
plot(f,10*log10(X*winGain/rms(win)),'LineWidth',1,'Color','b')
plot(f,10*log10(X),'LineWidth',1,'Color','m')
axis([95 105 -10 0])
xlabel('Frequency [Hz]')
ylabel('Power Spectrum [dB re 1V]')
legend({'W_{NG}','W_{CG}'},...
    'Location','NorthEast')
box on
set(gcf,'PaperPositionMode','auto')
title('Window Gain Effect (CG, Detail of 100 Hz tone)')
print('Window Gain Effect (CG, Detail)','-dpng','-r250')
savefig('Window Gain Effect (CG, Detail)')

% Window Gain Effect (Transient)
nfft = 8192;
delaySamples = round(nSamples/20);
winhan = hamming(nSamples);
winrec = rectwin(nSamples);
x_han1 = x_gun.*winhan;
x_han2 = circshift(x_gun,delaySamples).*winhan;
x_rec = x_gun.*winrec;
winGain1 = rms(winhan);
winGain2 = rms(x_han1)/rms(x_gun);
winGain3 = rms(x_han2)/rms(x_gun);
X_han1 = fft(x_han1,nfft)/winGain2; % original FFT
X_han1 = X_han1.*conj(X_han1)/(nSamples*fs);
X_han1 = [X_han1(1); 2*X_han1(2:ceil(nfft/2))];
X_han2 = fft(x_han2,nfft)/winGain3; % original FFT
X_han2 = X_han2.*conj(X_han2)/(nSamples*fs);
X_han2 = [X_han2(1); 2*X_han2(2:ceil(nfft/2))];
X_rec = fft(x_rec,nfft); % original FFT
X_rec = X_rec.*conj(X_rec)/(nSamples*fs);
X_rec = [X_rec(1); 2*X_rec(2:ceil(nfft/2))];
f = (0:length(X)-1)*fs/(2*length(X));

figure
hold on
plot(f,10*log10(X_rec),'LineWidth',1,'Color','b')
plot(f,10*log10(X_han1),'LineWidth',1,'Color','m')
plot(f,10*log10(X_han2),'LineWidth',1,'Color','g')
xlabel('Frequency [Hz]')
ylabel('Power Spectral Density [dB re 1V]')
legend({'Rectwin','Hanning (Delay = 0 ms)'...
    sprintf('Hanning (Delay = %0.0f ms)',delaySamples*1e3/fs)},...
    'Location','NorthEast')
box on
xlim([0 250])
set(gcf,'PaperPositionMode','auto')
title('Window Gain Effect (TG, Spectrum)')
print('Window Gain Effect (TG, Spectrum)','-dpng','-r250')
savefig('Window Gain Effect (TG, Spectrum)')

figure
hold on
plot(t,x_rec/max(x_rec),'LineWidth',1,'Color','b')
plot(t,x_han1/max(x_han1),'LineWidth',1,'Color','m')
plot(t,x_han2/max(x_han2),'LineWidth',1,'Color','g')
xlabel('Time [s]')
ylabel('Normalised Amplitude')
legend({'Rectwin','Hanning (Delay = 0 ms)'...
    sprintf('Hanning (Delay = %0.0f ms)',delaySamples*1e3/fs)},...
    'Location','NorthEast')
box on
ylim([-1.5 1.5])
set(gcf,'PaperPositionMode','auto')
title('Window Gain Effect (TG, Waveform)')
print('Window Gain Effect (TG, Waveform)','-dpng','-r250')
savefig('Window Gain Effect (TG, Waveform)')


