# unsure if this is correct

import numpy as np
import matplotlib.pyplot as plt

# Create signal
srate = 1000
time = np.arange(-3, 3.00001, 1/srate)
pnts = len(time)
freqmod = np.exp(-time**2)*10+10
freqmod = freqmod + np.linspace(0, 10, pnts)
signal = np.sin(2*np.pi * ((time + np.cumsum(freqmod))/srate))

# Plot the signal
fig, axs = plt.subplots(3, 1, figsize=(8, 6))

axs[0].plot(time, signal, linewidth=1)
axs[0].set_xlabel('Time (s)')
axs[0].set_title('Time-domain signal')

n = 500
hz = np.linspace(0, srate, n)
tf = np.zeros((int(pnts/n)-1, len(hz)))
tv = np.zeros(int(pnts/n)-1)

for i in range(int(pnts/n)-1):
    datasnip = signal[i*n:(i+1)*n]
    pw = np.abs(np.fft.fft(datasnip))**2 / n  # Normalization
    tf[i, :len(hz)] = pw[:len(hz)]
    tv[i] = np.mean(time[i*n:(i+1)*n])

axs[1].imshow(tf.T, aspect='auto', origin='lower', extent=[tv.min(), tv.max(), 0, 40], cmap='hot', vmin=0, vmax=np.max(tf)/5)
axs[1].set_xlabel('Time (s)')
axs[1].set_ylabel('Frequency (Hz)')
axs[1].set_title("Professor's Solution")

# My Solution
import scipy.io

data = scipy.io.loadmat('spectral\spectral_codeChallenge.mat')
signal = data['signal'].flatten()

n = len(signal)
winlength = int(srate / 2)
step = int(winlength / 2)  # 50% overlap
td_wlen_start = np.arange(0, n-winlength, step)  # Time start of each window

hzW = np.linspace(0, srate/2, int(winlength/2))
signalpowW = []

for wi in range(len(td_wlen_start)):
    datachunk = signal[td_wlen_start[wi]:td_wlen_start[wi]+winlength]
    tmppow = np.abs(np.fft.fft(datachunk)/winlength)**2
    signalpowW.append(tmppow[:int(winlength/2)])

signalpowW = np.array(signalpowW).T

# Create a time vector based on the window start times
time_vector = (td_wlen_start) / srate

axs[2].imshow(signalpowW, aspect='auto', origin='lower', extent=[time_vector.min(), time_vector.max(), 0, 40], cmap='viridis', vmin=0, vmax=np.max(signalpowW)/5)
axs[2].set_xlabel('Time (s)')
axs[2].set_ylabel('Frequency (Hz)')
axs[2].set_title('My Solution')
fig.colorbar(axs[2].images[0], ax=axs[2])

plt.tight_layout()
plt.show()