import wfdb
import ECG_BASELINE_FCNS
import matplotlib.pyplot as plt
import numpy as np

# number of samples to read
samples = 10000

fileName = np.array([100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 111, 112, 113, 114, 115, 116, 117, 118, 119, 121, 122, 123, 124, 200, 201, 202, 203, 205, 207, 208, 209, 210, 212, 213, 214, 215, 217, 219, 220, 221, 222, 223, 228, 230, 231, 232, 233, 234])

for i in range(0, len(fileName)):
    signals, fields = wfdb.srdsamp('in_samples/'+str(fileName[i]), sampto=samples)

    # get sampling rate
    fs = fields['fs']

    # create 1D numpy array of ECG data for filter usage
    ecgsignal = (np.transpose(signals))[0]

    # get signal time array (for display)
    tm = np.arange(0, len(ecgsignal)/fs, 1/fs)

    butterworthSignal = ECG_BASELINE_FCNS.ECG_BASELINE(ecgsignal, fs, 'butterworth')
    movAvSampleSignal = ECG_BASELINE_FCNS.ECG_BASELINE(ecgsignal, fs, 'movavSample')
    movAvTimeSignal = ECG_BASELINE_FCNS.ECG_BASELINE(ecgsignal, fs, 'movavSample')
    savgolSignal = ECG_BASELINE_FCNS.ECG_BASELINE(ecgsignal, fs, 'savitzkygolay')
    waveletSignal = ECG_BASELINE_FCNS.ECG_BASELINE(ecgsignal, fs, 'wavelet')

    np.save('out_butterworth/'+str(fileName[i])+'.npy', butterworthSignal)
    np.save('out_movAverage/' + str(fileName[i]) + '.npy', movAvSampleSignal)
    np.save('out_savGol/' + str(fileName[i]) + '.npy', savgolSignal)
    np.save('out_wavelet/' + str(fileName[i]) + '.npy', waveletSignal)

    np.save('out_signalInfo/' + str(fileName[i]) + 'freq.npy', fs)
    np.save('out_signalInfo/' + str(fileName[i]) + 'time.npy', tm)

    # plt.figure(1)
    # plt.subplot(511)
    # plt.plot(tm, ecgsignal)
    # plt.title('Original ECG signal')
    #
    # plt.subplot(512)
    # plt.plot(tm, butterworthSignal)
    # plt.title('Butterworth filter')
    #
    # plt.subplot(513)
    # plt.plot(tm, movAvTimeSignal)
    # plt.title('Moving average filter')
    #
    # plt.subplot(514)
    # plt.plot(tm, savgolSignal)
    # plt.title('Savitzky - Golay filter')
    #
    # plt.subplot(515)
    # plt.plot(tm, waveletSignal)
    # plt.title('Wavelet adaptive filter')
    #
    # plt.show()