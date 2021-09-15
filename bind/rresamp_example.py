#!/usr/bin/env python3
import sys
sys.path.extend(['.','..'])
import numpy as np
import liquid as dsp
import matplotlib.pyplot as plt

# options
P, Q, m = 3, 5, 20
num_blocks = 30

# generate pulse
num_samples = num_blocks * P
t = np.arange(2*m+1) - m
x = np.zeros((num_samples,), dtype=np.csingle)
x[:2*m+1] = np.sinc(0.2*t) * np.hamming(2*m+1)

# create resampler
resamp = dsp.rresamp(P, Q, m)
print(resamp)

# run resampler on pulse
y = resamp.execute(x)

# compute spectral responses
nfft = 2400
X = 20*np.log10(np.abs(np.fft.fftshift(np.fft.fft(x,   nfft))))
Y = 20*np.log10(np.abs(np.fft.fftshift(np.fft.fft(y/resamp.rate, nfft))))

# compute time and frequency arrays
tx = np.arange(len(x))
ty = np.arange(len(y))/resamp.rate - resamp.delay
fx = np.arange(nfft)/nfft-0.5
fy = fx * resamp.rate

# plot time and frequency responses
fig, (ax1, ax2) = plt.subplots(2,figsize=(8,8))
ax1.plot(tx, np.real(x), ty, np.real(y))
ax1.set_xlabel('Delay [samples]')
ax1.set_ylabel('Impulse Response')
ax1.grid(True, zorder=5)
ax1.legend(('input','output'))
ax2.plot(fx, X, fy, Y)
ax2.set_xlabel('Normalized Frequency [f/F_s]')
ax2.set_ylabel('Power Spectral Density [dB]')
ax2.set(xlim=(-0.5,0.5))
ax2.grid(True, zorder=5)
plt.show()

