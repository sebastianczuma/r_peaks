import numpy as np
import R_PEAKS
import ECG_BASELINE
import wfdb
import matplotlib.pyplot as plt

#which = 'butterworth'
which = 'movAverage'
#which = 'savGol'
#which = 'wavelet'

'-----------------------------------------------------------------------------'

#file_name = np.array([100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 111, 112, 113, 114, 115, 116, 117, 118, 119, 121, 122, 123, 124, 200, 201, 202, 203, 205, 207, 208, 209, 210, 212, 213, 214, 215, 217, 219, 220, 221, 222, 223, 228, 230, 231, 232, 233, 234])
##file_name = np.array([108])
#print(len(file_name))
##ecg_filtered = np.empty(0)
#for i in range(0, len(file_name)):
#	signals, fields = wfdb.srdsamp('Signals/' + str(file_name[i]))
#	ecgsignal = (np.transpose(signals))[0]
#	fs = fields['fs']
#	r_val, r_index = R_PEAKS.FUNC_r_detection(ecgsignal, fs)
#	print(str(file_name[i]) + ' ' + str(len(r_index)))
#	print(len(r_index))
##   ecg_filtered_p = np.load('ECG_BASELINE/ref/out_'+str(which)+'/'+str(file_name[i])+'.npy')
##   ecg_filtered = np.append(ecg_filtered, ecg_filtered_p)
##   fs = np.load('ECG_BASELINE/ref/out_signalInfo/'+str(file_name[i])+'freq.npy')
    
    #t = np.load('ECG_BASELINE/ref/out_signalInfo/'+str(file_name[i])+'time.npy')
    #milivolty
    #ecg_filtered = ecg_filtered*1000

    #r_value, r_index = R_PEAKS.FUNC_r_detection(ecg_filtered, fs)

    #np.save('R_PEAKS/ref/R_out_'+str(which)+'/R_index/r_index'+str(file_name[i])+'.npy', r_index)
    #np.save('R_PEAKS/ref/R_out_'+str(which)+'/R_value/r_value'+str(file_name[i])+'.npy', r_value)

#r_index = np.load('R_PEAKS/ref/R_out_'+str(which)+'/R_index/r_index'+str(file_name[i])+'.npy')
#r_value = np.load('R_PEAKS/ref/R_out_'+str(which)+'/R_value/r_value'+str(file_name[i])+'.npy')

#R_PEAKS.PRINT_r(ecg_filtered, r_index, r_value)

#samples=10000;
# 1000 - 11000 , 222
# 1000 - 3000
# 2000 - 3000
signals, fields = wfdb.srdsamp('Signals/222')#, sampfrom=0,sampto=11000)

# get sampling rate
fs = fields['fs']

# create 1D numpy array of ECG data for filter usage
ecgsignal = (np.transpose(signals))[0]

#ecgb = R_PEAKS.FUNC_flat_izo(ecgsignal)

#ecgb = ecgb * ecgb * ecgb * ecgb * 100

r_values, r_index = R_PEAKS.FUNC_r_detection(ecgsignal, fs)

#r_values = np.empty(len(r_index))

#for i in range(len(r_index)):
#	z = int(r_index[i])
#	r_values[i] = ecgsignal[z]


R_PEAKS.PRINT_r(ecgsignal, r_index, r_values)
#R_PEAKS.PRINT_any(ecgb)

##R_PEAKS.PRINT_all(ecg_filtered, fs)

#R_PEAKS.PRINT_all(ecgsignal, fs)


'---Uzycie funkcji FUNC_r_detection---'

#print(len(r_index))
#R_PEAKS.PRINT_r(ecgsignal, r_index, r_val)


# record = wfdb.rdsamp("Signals/222", channels=[0])
# d_signal = record.adc()[:,0]
# #peak_indices = wfdb.processing.gqrs_detect(x=d_signal, fs=record.fs, adcgain=record.adcgain[0], adczero=record.adczero[0], threshold=1.0)
# min_bpm = 10
# max_bpm = 350
# min_gap = record.fs*60/min_bpm
# max_gap = record.fs*60/max_bpm
# #new_indices = wfdb.processing.correct_peaks(d_signal, peak_indices=peak_indices, min_gap=min_gap, max_gap=max_gap, smooth_window=150)

# r_val, r_index = R_PEAKS.FUNC_r_detection(d_signal, record.fs)

# r_inx = [len(r_index)]
# for i in range(len(r_index)):
# 	r_inx.append(int(r_index[i]))

# new_indices_nn = wfdb.processing.correct_peaks(d_signal, peak_indices=r_inx, min_gap=min_gap, max_gap=max_gap, smooth_window=150)

# #print(new_indices)
# plt.figure(5)
# plt.plot(d_signal)
# #plt.plot(new_indices, d_signal[new_indices], marker='x', color='r', ls='')
# plt.plot(new_indices_nn, d_signal[new_indices_nn], marker='o', color='y', ls='')
# plt.xlabel('zalamki R')

# plt.show()

# print(len(new_indices_nn))


#signal = R_PEAKS.FUNC_filters(ecgsignal, fs)
#R_PEAKS.PRINT_any(R_PEAKS.FUNC_running_mean(R_PEAKS.FUNC_flat_izo(ecgsignal)))

#signal = R_PEAKS.FUNC_flat_izo(ecgsignal)
#ecg_b = R_PEAKS.FUNC_filter(signal, fs)
#diff = R_PEAKS.FUNC_diff(ecg_b)
#sqr = R_PEAKS.FUNC_sqr(diff)
#ecgmf = R_PEAKS.FUNC_signal_integration(sqr, fs)
#left_T, right_T = R_PEAKS.FUNC_prepare_for_maxima_search(ecgmf, fs)
#r_value, r_index = R_PEAKS.FUNC_find_max(left_T, right_T, ecgsignal)

#R_PEAKS.PRINT_compare(ecgsignal, ecg_b)


