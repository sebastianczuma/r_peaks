import matplotlib.pyplot as plt
import scipy.signal as signal
import numpy as np
import pywt as pywt #leave it here even if not yet used
from statsmodels.robust import mad

# brief: gui-like function to call different filters
# function's name: ECG_BASELINE
# parameters:
#               ecgSignal   - numpy 1D array with ECG data
#               fs          - sampling rate
#               method      -   'butterworth': butterworth cascade filter
#                               'movavSample': moving average filter with number of samples length window
#                               'movavTime': moving average filter with time length window
#                               'savitzkygolay': Savitzky - Golay filter
# return value:
#               filteredSignal - numpy 1D array with denoised ECG data
#
# description:  function calculate denoised ECG signal with chosen method and author-proposed parameters
def ECG_BASELINE(ecgSignal, fs, method):

    if method == 'butterworth':

        filteredSignal = BUTTERWORTH(ecgSignal, fs)

    elif method == 'movavSample':

        filteredSignal = MOVINGAVERAGE(ecgSignal, 100, fs, 'sample')

    elif method == 'movavTime':

        filteredSignal = MOVINGAVERAGE(ecgSignal, 0.1, fs, 'time')

    elif method == 'savitzkygolay':

        filteredSignal = SAVITZKYGOLAY(ecgSignal, 3, 9)

    elif method == 'wavelet':

        filteredSignal = WAVELET(ecgSignal)

    else:

        filteredSignal = 0

    return filteredSignal

# brief: calculate signal filtered by Butterwoerth filter
# function's name: BUTTERWORTH
# parameters:
#               ecgSignal   - numpy 1D array with ECG data
#               fs          - sampling rate
#               response    - 'true': plot filter response, 'false': do not plot signal response (default)
# return value:
#               filteredSignal - numpy 1D array with denoised ECG data
#
# description:  function uses two cascade butterworth filters to eliminate ECG noises. First lowpass filter to
#               cut frequencies over 40Hz, next, highpass filter to cut frequencies below 0.5Hz
def BUTTERWORTH(ecgSignal, fs):

    order = 2

    # hum noise filter
    f1 = 59
    f2 = 61
    b1, a1 = signal.butter(order, [f1/(fs/2), f2/(fs/2)], btype='bandstop')
    filteredSignal = signal.filtfilt(b1, a1, ecgSignal)

    # min ecg freq
    f3 = 0.8
    b2, a2 = signal.butter(order, f3/(fs/2), btype='highpass')
    filteredSignal = signal.filtfilt(b2, a2, filteredSignal)

    # muscle noise filter
    f4 = 35;
    b3, a3 = signal.butter(order, f4/(fs/2), btype='lowpass')
    filteredSignal = signal.filtfilt(b3, a3, filteredSignal)

    return filteredSignal

# brief: calculate signal filtered by moving average filter
# function's name: MOVINGAVERAGE
# parameters:
#               ecgSignal   - numpy 1D array with ECG data
#               window      - number of samples in window or window time length in [s]
#               fs          - sampling rate
#               method      - 'time': window in time domain; 'sample': window in number of samples domain (default)
# return value:
#               filteredSignal - numpy 1D array with denoised ECG data
#
# description:  function uses moving average filter to smooth ECG signal
def MOVINGAVERAGE(ecgSignal, window, fs, method='sample'):

    movAvSignal = np.zeros(len(ecgSignal))

    if method == 'time':
        print(np.around((window * fs)/2))
        window = int(np.around((window * fs)/2))
    else:
        window = int(np.around(window/2))

    for s in range (1, len(ecgSignal)):

        if((s <= window) or (len(ecgSignal) - s <= window)):
            movAvSignal[s-1] = ecgSignal[s-1]
        else:
            movAvSignal[s-1] = (1/(2 * window + 1))*np.sum(ecgSignal[(s - window - 1):(s + window - 1)])

        filteredSignal = np.transpose(ecgSignal - movAvSignal)

    return filteredSignal

# brief: calculate signal filtered by Savitzky - Golay filter
# function's name: SAVITZKYGOLAY
# parameters:
#               ecgSignal   - numpy 1D array with ECG data
#               order       - order of filter
#               window      - window length
# return value:
#               filteredSignal - numpy 1D array with denoised ECG data
#
# description:  function uses Savitzky - Golay filter to smooth ECG signal
def SAVITZKYGOLAY(ecgSignal, order, window):

    filteredSignal = signal.savgol_filter(ecgSignal, window, order)
    return filteredSignal

# brief: calculate signal filtered by using wavelet filter
# function's name: WAVELET
# parameters:
#               ecgSignal   - numpy 1D array with ECG data
# return value:
#               filteredSignal - numpy 1D array with denoised ECG data
#
# description:  function uses wavelet transformation to filter signal.
#               First, wavelet transformation is calculated to get
#               wavelet base of signal. Than, using hard thresholding,
#               vector base is 'reduced' and signal is reconstructed.
def WAVELET(ecgSignal):

    coeffs = pywt.wavedec(ecgSignal, 'db5', level=6)

    sigma = mad(coeffs[-1])
    uthresh = sigma * np.sqrt(2 * np.log(len(coeffs)))
    denoised = coeffs[:]

    denoised[1:] = (pywt.threshold(i, value=uthresh, mode='hard') for i in denoised[1:])

    filteredSignal = pywt.waverec(denoised, 'db5')
    return filteredSignal
