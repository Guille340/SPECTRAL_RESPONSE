% [X,f] = spectralResponse(x,fs,varargin)
%
% DESCRIPTION
% Calculates the narrowband spectral response of input waveform X using a
% specific window and number of DFT bins. The spectrum or the spectral
% density for the metrics Energy(RMS^2*LENGTH(X)), Power (RMS^2), and Exposure
% (RMS^2*LENGTH(X)/FS).
%
% INPUT VARIABLES
% - x: input waveform.
% - fs: sampling rate [Hz]
% - PROPERTY/VALUE (varargin(n,n+1)): the following propeties are
%   available. In the function call, type the property string followed
%   by its value (comma-separated). Property/value pairs are variable
%   input arguments, and must be specified after the fixed arguments.
%   One or more different property/value pairs can be specified in the 
%   function call. 
%   ¬ 'NumberOfBins': number of DFT bins. For efficiency, use a power of 2.
%      If the number is higher than LENGTH(X), the function will apply zero-
%      padding to the input signal. If the number is lower than LENGTH(X), the
%      input signal will be wrapped using DATAWRAP. By DEFAULT, the smallest
%      power of 2 larger than LENGTH(X) is used.
%   ¬ 'Window': numeric vector representing the weighting window in time.
%      It must be the same length as input X. See WINDOW function. By
%      DEFAULT a rectangular window ONES(LENGTH(X)) is used (see RECTWIN).
%   ¬ 'Spectrum': string indicating the type of spectral response to be 
%      processed. There are six options available (Note: G = window gain, 
%      u = units of input X):
%      ~ 'pow': power autospectrum. FFT(X).^2./(G^2*LENGTH(X)^2) [u^2]
%      ~ 'ene': energy autospectrum. FFT(X).^2./(G^2*LENGTH(X)) [u^2]
%      ~ 'exp': exposure autospectrum. FFT(X).^2./(G^2*LENGTH(X)*FS) [u^2*s]
%      ~ 'psd': power spectral density. FFT(X).^2./(G^2*LENGTH(X)*FS) [u^2/Hz]
%         This is the DEFAULT option.
%      ~ 'esd': energy spectral density. FFT(X).^2./(G^2*FS) [u^2/Hz]
%      ~ 'xsd': exposure spectral density. FFT(X).^2./(G^2*FS^2) [u^2*s/Hz]
%   ¬ 'Signal': string indicating the type of time signal. This option will
%      affect how the window gain is calculated. Three options:
%      ~ 'random': stationary non-deterministic signal (e.g. white noise).
%         This is the DEFAULT option.
%      ~ 'coherent': stationary deterministic signal (e.g. tonal noise)
%      ~ 'transient': non-stationary pulsed or impulsive signal (e.g.
%         piling). Its spectral characteristics will vary with time.
%   ¬ 'FrequencyRange': string indicating the frequency region over which
%      the spectrum is calculated. There are three options available:
%      ~ 'onesided': single-sided spectrum (0 to FS/2). The power of all
%         frequency samples is doubled, except for the first one (DC). This
%         is the DEFAULT option.
%      ~ 'twosided': double-sided spectrum (0 to FS).
%      ~ 'centered': double-sided spectrum (-FS/2 to FS/2).
%   
% OUTPUT VARIABLES
% - X: narrowband power spectrum processed according to the input arguments.
% - f: vector of frequencies for X [Hz]
%
% CONSIDERATIONS & LIMITATIONS
% - SPECTRALRESPONSE produces the same results as PERIODOGRAM, except that
%   the former has a number of useful advantages. In particular:
%   ¬ The exposure metric can be represented (not only the power and
%   energy).
%   ¬ The type of window gain correction can be controlled by selecting
%     the signal type ('random','coherent','transient'). PERIODOGRAM uses
%     by default a coherent gain for 'power' and 'energy' options and
%     a noise-like gain for 'psd' and 'esd' options.
%   ¬ The code is implemented in a simple way for it to help understand
%     how the different narrowband spectra are calculated.
%
% FUNCTION DEPENDENCIES
% - None
%
% FUNCTION CALL
% 1) [X,f] = spectralResponse(x,fs)
% 2) [X,f] = spectralResponse(...,PROPERTY,VALUE)
%     
%    PROPERTY              VALUE
%    --------------------------------------------
%    'NumberOfBins'        Number of Fft bins
%    'Window'              Vector same length as x
%    'Spectrum'            'power','energy','psd','esd'
%    'Signal'              'random','coherent','transient'
%    'FrequencyRange'      'onesided','twosided','centered'

% VERSION 1.1
% Date: 22 May 2021
% Revised: Guillermo Jimenez Arranz
% - Updated the help.
%
% VERSION 1.0
% Date: 04 May 2021
% Author: Guillermo Jimenez Arranz
% Email: gjarranz@gmail.com

function [X,f] = spectralResponse(x,fs,varargin)

% INPUT ARGUMENTS
% Verify number of Input Arguments
narginchk(2,10)
if rem(nargin-2,2)
    error('Variable input arguments must come in pairs (PROPERTY,VALUE)')
end
% Verify Compulsory Input Parameters
if ~isnumeric(x)
    error('Input argument X must be a numeric vector or array')
end
if ~isnumeric(fs) || fs < 0
    fs = 2;
    warning(['Input argument FS must be a numeric value higher than 0. '...
        'FS = 2 will be used'])
end
% Parameters
x = x(:);
nSamples = length(x);
% Extract and Verify Input Properties
validProperties = {'numberofbins','window','spectrum',...
    'signal','frequencyrange'};
properties = lower(varargin(1:2:end));
if any(~ismember(properties,validProperties))
    error('One or more PROPERTY is not recognised')
end
% Default Input Values
nfft = 2^ceil(log2(nSamples));
win = ones(nSamples,1);
specType = 'psd';
signalType = 'random';
freqRange = 'onesided';
% Extract and Verify Input Values
values = varargin(2:2:end);
nPairs = (nargin - 1)/2; % number of (PROPERTY,VALUE) pairs
for m = 1:nPairs
    property = properties{m};
    switch property % populate with more properties if needed
        case 'numberofbins'
            nfft = values{m};
            if ~isnumeric(nfft) || numel(nfft) > 1 || nfft < 1
                nfft = 2^ceil(log2(nSamples));
                warning(['Non-supported value for PROPERTY = '...
                    '''NumberOfBins''. ''NumberOfBins'' = %d '...
                    'will be used'],nfft)
            end
        case 'window'
            win = values{m};
            if ~isnumeric(win) || ~isvector(win) || length(win) ~= nSamples
                win = ones(nSamples,1);
                warning(['Input property ''Window'' must be a numeric '...
                    'vector the same length as X. ''Window'' = '...
                    'ONES(LENGTH(X),1) will be used'])
            end   
        case 'spectrum'
            specType = lower(values{m});
            if ~ismember(specType,{'pow','ene','exp','psd','esd','xsd'})
                specType = 'psd';
                warning(['Non-supported value for PROPERTY = '...
                    '''Spectrum''. ''Spectrum'' = ''psd'' '...
                    'will be used'])
            end  
        case 'signal'
            signalType = lower(values{m});
            if ~ismember(signalType,{'random','coherent','transient'})
                signalType = 'random';
                warning(['Non-supported value for PROPERTY = '...
                    '''Signal''. ''Signal'' = ''%s'''...
                    'will be used'],signalType)
            end
        case 'frequencyrange'
            freqRange = lower(values{m});
            if ~ismember(freqRange,{'onesided','twosided','centered'})
                freqRange = 'onesided';
                warning(['Non-supported value for PROPERTY = '...
                    '''FrequencyRange''. ''FrequencyRange'' = '...
                    '''onesided''' 'will be used'])
            end
    end 
end

% Window Input Signal
if isvector(x), x = x(:); end % if x is a vector, make it a column
xwin = x.*win;

% Calculate Window Gain
switch signalType
    case 'random'
        winGain = rms(win);        
    case 'coherent'
        winGain = mean(win);
    case 'transient'
        winGain = rms(xwin)/rms(x);  
end

% Calculate Normalising Factor
switch specType
    case 'ene'
        normFactor = nSamples;
    case 'pow'
        normFactor = nSamples*nSamples;
    case 'exp'
        normFactor = nSamples*fs;
    case 'esd'
        normFactor = fs;
    case 'psd'                     
        normFactor = nSamples*fs;
    case 'xsd'
        normFactor = fs*fs;
end

% Normalised Spectrum
xwin = datawrap(xwin,nfft);
xfft = fft(xwin)/winGain;
switch freqRange
    case {'twosided','centered'}
        X = xfft.*conj(xfft)/normFactor;
        f = (0:fs/nfft:fs-fs/nfft);
        if strcmp(freqRange,'centered')
            X = fftshift(X); 
            f = f - fs/2 + fs/nfft;
        end        
    case 'onesided'
        halfPoint = ceil((nfft + 1)/2);
        xfft = xfft(1:halfPoint);
        X = xfft.*conj(xfft)/normFactor;
        X(2:end,:) = 2*X(2:end,:); 
        f = (0:fs/nfft:fs/2);   
end
