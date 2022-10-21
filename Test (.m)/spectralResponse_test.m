% Input
sourceName = 'AKP_3400LLX_080_2000_80';
sourceNameTitle = 'AKP 3400LLX 080 2000 80';
fileName = 'AKP_3400LLX_080_2000_80 (fs = 2kHz).csv';

% Extract Signal
csvData = csvread(fileName,1,0);
t = csvData(:,1)*1e-3; % time [s]
x = csvData(:,2); % pressure [bar m]
sampleInterval = unique(round(diff(t)*1e6)*1e-6);

% General Parameters
pref = 1e-6;
bar2pascal = 1e5;
fs = 1/sampleInterval;

% % Signal
% fs = 2000;
% nSamples = 4001;
% t = (0:nSamples - 1)/fs;
% f = [100.09765625 250];
% x = sin(2*pi*f(1)*t) + 0.5*sin(2*pi*f(2)*t) ;

x = x + 0.5;

% % Filter Signal
[X_psd,f] = spectralResponse(x,fs,'NumberOfBins',8192,'Spectrum','psd',...
    'Window',hamming(length(x)),'FrequencyRange','centered');
[Xp_psd,fp] = periodogram(x,hamming(length(x)),8192,fs,'psd','centered');

% PERIODOGRAM PSD
% MATLAB's periodogram 'psd' uses NG = rms(win)
% to normalise the spectrum by the window. The two lines above
% produce identical results for any window. That NG dividing factor
% does not properly compensate for the effect of the window when
% the signal being windowed is transitory, and only works with
% random continuous noise. The only way of fixing that is using a
% window gain for transitory signals, calculated as TG =
% rms(x.*win)/rms(x). This gain gives much better results than the
% gain for coherent signals CG = mean(win).

[X_pow,f] = spectralResponse(x,fs,'NumberOfBins',8192,'Spectrum','pow',...
    'Window',blackman(length(x)));
[Xp_pow,fp] = periodogram(x,hamming(length(x)),8192,fs,'power');

% PERIODOGRAM POWER
% MATLAB's periodogram 'power' uses CG = mean(win) to normalise
% the spectrum by the window. It also multiplies by corrFactor^2,
% rather than just corrFactor.


