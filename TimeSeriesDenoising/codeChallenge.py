import numpy as np
import matplotlib.pyplot as plt

# Load the signal package
# (No need for this in Python)

N = 4000

# Generate the original signal
t = np.linspace(0, 10 * np.pi, N)
origSignal = np.sin(t) * np.linspace(-1, 1, N) + np.random.randn(N)
r = np.random.permutation(N)
nn = int(np.round(N * 0.05))
origSignal[r[:nn]] = (1 + np.random.rand(nn)) * 10
origSignal[r[-nn:]] = -(1 + np.random.rand(nn)) * 10

# Initialize cleanedSignal with origSignal
cleanedSignal = origSignal.copy()

# Remove positive noise spikes
p2r = np.where(origSignal > 5)[0]
k = 5
for i in p2r:
    start = max(0, i - k)
    end = min(N, i + k + 1)
    cleanedSignal[i] = np.median(cleanedSignal[start:end])

# Remove negative noise spikes
p2r = np.where(origSignal < -5)[0]
k = 5
for i in p2r:
    start = max(0, i - k)
    end = min(N, i + k + 1)
    cleanedSignal[i] = np.median(cleanedSignal[start:end])

# Mean-smooth
k = 150
for i in range(N):
    start = max(0, i - k)
    end = min(N, i + k + 1)
    cleanedSignal[i] = np.mean(cleanedSignal[start:end])

# Plot the original signal and the lecturer's solution
plt.figure()
plt.subplot(3, 1, 1)
plt.plot(origSignal, linewidth=2)
plt.title("Original Signal")

plt.subplot(3, 1, 2)
plt.plot(cleanedSignal, linewidth=2)
plt.title("Lecturer's Solution")

# Your Solution
outlierThreshold = 5
k = 5
myCleanedSignal = origSignal.copy()

# Positive outliers
outlierIndices = np.where(origSignal > outlierThreshold)[0]
for i in outlierIndices:
    start = max(0, i - k)
    end = min(N, i + k + 1)
    myCleanedSignal[i] = np.median(origSignal[start:end])

# Negative outliers
outlierIndices = np.where(origSignal < -outlierThreshold)[0]
for i in outlierIndices:
    start = max(0, i - k)
    end = min(N, i + k + 1)
    myCleanedSignal[i] = np.median(origSignal[start:end])

# Moving average smoothing
windowSize = 151
movingAvgFilter = np.ones(windowSize) / windowSize
mySmoothSignal = np.convolve(myCleanedSignal, movingAvgFilter, mode='same')

plt.subplot(3, 1, 3)
plt.plot(mySmoothSignal, color=[0.8, 0.2, 0.1], linewidth=2)
plt.title('Your Solution')

plt.tight_layout()
plt.show()