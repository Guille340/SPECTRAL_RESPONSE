Spectral Response
=======================

MATLAB code for calculating the spectral response of a discrete time signal.
The spectral response can be calculated for three different metrics (power,
energy, exposure) and as two different representations (energy or density).
Therefore, six types of spectra can be obtained: PS (power spectrum),
ES (energy spectrum), XS (exposure spectrum), PSD (power spectral density),
ESD (energy spectral density), and XSD (exposure spectral density).

The function SPECTRALRESPONSE has similar functionality to MATLAB's PERIODOGRAM,
except that with SPECTRALRESPONSE the type of window weighting can be selected
(Noise Gain NG or Coherent Gain CG), and includes the exposure spectrum as an
extra option (particularly useful for underwater acoustics, where sound exposure
is a common metric). Furthermore, SPECTRALRESPONSE has been written with an
educational purpose, to learn how to properly calculate the spectral response
and understand the processing steps.

For further details on FFT analysis and spectrum calculations check the "Docs"
folder. The "how-to-use-the-fft-and-matlab-s-pwelch-function.pdf" article
is particularly useful. I made a summary of the different ways of calculating
a signal's narrowband in "Notes on Power Spectrum".

[Guillermo Jiménez Arranz, 16 Jun 2021]





