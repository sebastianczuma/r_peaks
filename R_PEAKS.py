import matplotlib.pyplot as plt
import scipy.signal as ss
import numpy as np
import math

window_length = 0.15 #sec
'------------------------------------------------------------------------------'


# opis: wyrownuje izolinie.
# nazwa funkcji: FUNC_flat_izo
# parametry:
#               ecg_filtered - numpy 1D array
# zwraca:
#               izo - numpy 1D array
def FUNC_flat_izo(ecg_filtered):
	# mV
	ecg_filtered = ecg_filtered * 1000

	N = len(ecg_filtered)

	delta_m = 3;
	delta_d = 3;

	max_g = np.empty(N)
	min_g = np.empty(N)
	h = np.empty(N)
	a = np.empty(N)
	n = np.empty(N)

	for i in range(0, N):
		if i==0:
			max_g[i] = ecg_filtered[i]
			min_g[i] = ecg_filtered[i]
		else:
			if ecg_filtered[i] > max_g[i-1]:
				max_g[i] = max_g[i-1] + delta_m * delta_d
			else:
				max_g[i] = max_g[i-1] - delta_d

			if ecg_filtered[i] < min_g[i-1]:
				min_g[i] = min_g[i-1] - delta_m * delta_d
			else:
				min_g[i] = min_g[i-1] + delta_d
		# V
		h[i] = (ecg_filtered[i] - (max_g[i] + min_g[i]) / 2) / 1000
		
	return h


# opis: srednia.
# nazwa funkcji: FUNC_running_mean
# parametry:
#               signal - numpy 1D array
#               fs - integer
# zwraca:
#               numpy 1D array
def FUNC_running_mean(signal, fs):
	N = 2 * fs
	return np.convolve(signal, np.ones((N,))/N)[(N-1):]


# opis: filtr gornoprzepustowy.
# nazwa funkcji: FUNC_filter
# parametry:
#               izo - numpy 1D array
#				fs - integer
# zwraca:
#               new_filt - numpy 1D array
def FUNC_filter(izo, fs):
	izo = izo - np.mean(izo)

	Wn = [5 * 2 / fs, 15 * 2 / fs]
	#Wn = 15 * 2 / fs
	N = 3
	[a,b] = ss.butter(N,Wn, 'bandpass')
	ecg_b = ss.filtfilt(a, b, izo)

	#new_filt = np.empty(N)
	#new_filt = new_filt[31:N-1] - filt[32:N]/32 + filt[16:N-16] - filt[15:N-17] + filt[1:N-31]/32

	return ecg_b


# opis: oblicza i zwraca pochodna sygnalu.
# nazwa funkcji: FUNC_diff
# parametry:
#               ecg_b - numpy 1D array
# zwraca:
#               diff - numpy 1D array
def FUNC_diff(ecg_b):
	N = len(ecg_b)
	'''POCHODNA---------------------------------------------------------------'''
	#diff = np.diff(ecg_filtered)
	diff = (2*ecg_b[5:N]+ecg_b[4:N-1]-ecg_b[2:N-3]-2*ecg_b[1:N-4])/8
	'''Wektor wynikowy krotszy o 5 elementow'''

	return diff


# opis: oblicza i zwraca potege sygnalu.
# nazwa funkcji: FUNC_sqr
# parametry:
#               diff - numpy 1D array 
# zwraca:
#               sqr - numpy 1D array
def FUNC_sqr(diff):
	'''POTEGA-----------------------------------------------------------------'''
	sqr = diff**2

	return sqr


# opis: calkuje sygnal w ruchomym oknie i zwraca go.
# nazwa funkcji: FUNC_signal_integration
# parametry:
#               sqr - numpy 1D array 
#				fs - integer
# zwraca:
#               ecgmf - numpy 1D array
def FUNC_signal_integration(sqr, fs):
	dt = 1/float(fs)
	'''CALKOWANIE W RUCHOMYM OKNIE--------------------------------------------'''
	'''Utworzenie okna'''

	window_size = window_length * fs
	window_size = math.ceil(window_size)

	window = np.ones(window_size)
	'''Calkowanie'''
	temp = ss.lfilter(window,1,sqr)
	ecgmf = ss.medfilt(temp,9)
	ecgmf = ecgmf * dt
	'''Usuniecie opoznienia filtru'''
	delay = math.ceil(len(window) / 2)
	ecgmf = ecgmf[delay:len(ecgmf)]

	return ecgmf


# opis: wyszukuje i zwraca przyblizone granice QRS
#		potrzebne do wyszukania zalamka R (oraz zalamka S).
# nazwa funkcji: FUNC_prepare_for_maxima_search
# parametry:
#               ecgmf - numpy 1D array
#               fs - integer
# zwraca:
#               left_T - numpy 1D array
#				right_T - numpy 1D array
def FUNC_prepare_for_maxima_search(ecgmf, fs):

	r_mean = FUNC_running_mean(ecgmf, fs)

	region_poz = np.empty(len(ecgmf))

	for i in range(len(ecgmf)):
		if ecgmf[i] > 1.5 * r_mean[i]:
			region_poz[i] = 1
		else:
			region_poz[i] = 0

	region_poz = np.array([region_poz])

	'''Uzupelnienie zerem'''
	region_poz_LR = np.insert(region_poz, 0, 0)
	region_poz_RL = np.append(region_poz, 0)

	'''SZUKANE MAKSIMOW-------------------------------------------------------'''
	deltaLR = np.diff(region_poz_LR)
	deltaRL = np.diff(region_poz_RL)

	'''Wyznaczenie granic segmentow'''
	left  = np.where(deltaLR == 1)
	right = np.where(deltaRL == -1)

	left_T = np.transpose(left)
	right_T = np.transpose(right)

	return left_T, right_T


# opis: znajduje i zwraca wartosci QRS na podstawie przyblizonych granic QRS.
# nazwa funkcji: FUNC_find_qrs_values
# parametry:
#               left_T - numpy 1D array
#				right_T - numpy 1D array
#				ecg_filtered - numpy 1D array
#				ecgmf - numpy 1D array
# zwraca:
#               qrs_left_values - numpy 1D array
#               qrs_right_values - numpy 1D array
#               qrs_left_values_ecgmf - numpy 1D array
#               qrs_right_values_ecgmf - numpy 1D array
def FUNC_find_qrs_values(left_T, right_T, ecg_filtered, ecgmf):
	qrs_left_values = np.empty(len(left_T))
	qrs_right_values = np.empty(len(left_T))
	qrs_left_values_ecgmf = np.empty(len(left_T))
	qrs_right_values_ecgmf = np.empty(len(left_T))

	ecg_filtered[12:len(ecg_filtered)]

	for i in range(0,len(left_T)):
		qrs_left_values[i] = ecg_filtered[left_T[i]]
		qrs_right_values[i] = ecg_filtered[right_T[i]]
		qrs_left_values_ecgmf[i] = ecgmf[left_T[i]]
		qrs_right_values_ecgmf[i] = ecgmf[right_T[i]]

	return qrs_left_values, qrs_right_values, qrs_left_values_ecgmf, qrs_right_values_ecgmf


# opis: znajduje lokalne maksima w przyblizonych granicach QRS.
#		Zwraca wartosci i indeksy wyszukanych maksimow.
# nazwa funkcji: FUNC_find_max
# parametry:
#               left_T - numpy 1D array
#				right_T - numpy 1D array
#				ecg_filtered - numpy 1D array
# zwraca:
#               max_value - numpy 1D array
#               max_index - numpy 1D array
def FUNC_find_max(left_T, right_T, ecg_filtered):
	max_index = np.empty(len(left_T))
	max_value = np.empty(len(left_T))
	#obciecie poczatku sygnalu oryginlanego w celu dopasownania indeksow po calkowaniu i pochodnej wzgledem oryginalnego
	#ecg_filtered = ecg_filtered[18:len(ecg_filtered)]

	#sygnal w trakcie calkownia jest scinany na koncu

	for i in range(0,len(left_T)):
		start = int(left_T[i])
		end = int(right_T[i])

		max_value[i] = ecg_filtered[start]
		max_index[i] = start

		for j in range(start,end):
			if ecg_filtered[j] > max_value[i]:
				max_value[i] = ecg_filtered[j]
				max_index[i] = j

		#max_index[i] = np.argmax(ecg_filtered[left_T[i]:right_T[i]])
		#max_index[i] = max_index[i]+left_T[i]
		#max_value[i] = ecg_filtered[max_index[i]]

	'''for i in range(0,len(left_T)):
		max_index[i] = np.argmax(ecg_filtered[left_T[i]:right_T[i]])
		max_index[i] = max_index[i]+left_T[i]
		max_value[i] = ecg_filtered[max_index[i]]'''

	return max_value, max_index


# opis: znajduje lokalne minima w przyblizonych granicach QRS.
#		Zwraca wartosci i indeksy wyszukanych minimow.
# nazwa funkcji: FUNC_find_min
# parametry:
#               r_index - numpy 1D array
#				right_T - numpy 1D array
#				ecg_filtered - numpy 1D array
# zwraca:
#               min_value - numpy 1D array
#               min_index - numpy 1D array
def FUNC_find_min(left_T, right_T, ecg_filtered):
	min_index = np.empty(len(left_T))
	min_value = np.empty(len(left_T))

	for i in range(0,len(left_T)):
		start = int(left_T[i])
		end = int(right_T[i])

		min_value[i] = ecg_filtered[start]
		min_index[i] = start

		for j in range(start,end):
			if ecg_filtered[j] < min_value[i]:
				min_value[i] = ecg_filtered[j]
				min_index[i] = j

	'''	

	for i in range(0,len(left_T)):
		min_index[i] = np.argmin(ecg_filtered[left_T[i]:right_T[i]])
		min_index[i] = min_index[i]+left_T[i]
		min_value[i] = ecg_filtered[min_index[i]]'''

	return min_value, min_index


# opis: sprawdza czasy pomiedzy kolejnymi zalamkami R. Jesli interwal
#       jest mniejszy niz 250 ms zalamek jest uznawany za niepoprawny.
# nazwa funkcji: FUNC_check_rr_intervals
# parametry:
#               r_index - numpy 1D array
#				r_value - numpy 1D array
#				fs - integer
# zwraca:
#               final_r_value - numpy 1D array
#               final_r_index - numpy 1D array
def FUNC_check_rr_intervals(r_index, r_value, fs, ecg_filtered):
	r_peaks_to_delete_compute = np.empty(len(r_index))
	del_index = 0
	range_ms = np.ceil(0.36 * fs)

	for i in range(len(r_index)-1):
		rr_interval = r_index[i] - r_index[i - 1]

		if rr_interval < range_ms:
			if r_value[i] < 0.5 * r_value[i - 1]:
				r_peaks_to_delete_compute[del_index] = i
				del_index = del_index + 1

	to_delete = np.empty(del_index)

	for i in range(del_index):
		to_delete[i] = int(r_peaks_to_delete_compute[i])

	final_r_index = np.delete(r_index, to_delete)
	final_r_value = np.delete(r_value, to_delete)

	rr_i = np.empty(len(r_index))

	for i in range(len(r_index)):
		if i > 1:
			rr_i[i] = r_index[i] - r_index[i-1]
			#print(rr_i[i])

			if i > 9:
				avg_rr = (rr_i[i-9]+rr_i[i-8]+rr_i[i-7]+rr_i[i-6]+rr_i[i-5]+rr_i[i-4]+rr_i[i-3]+rr_i[i-2]+rr_i[i-1])/8
				#print(avg_rr)
				if rr_i[i] > 1.66 * avg_rr:
					start = int(r_index[i-1]) + int(0.2 * fs)
					end = int(r_index[i]) - int(0.2 * fs)
					new_r_value = max(ecg_filtered[start:end])
					new_r_index = np.argmax(ecg_filtered[start:end]) + start
					#print(avg_rr)
					#print(end)
					if i < len(final_r_index):
						final_r_index = np.insert(final_r_index, i, new_r_index)
						final_r_value = np.insert(final_r_value, i, new_r_value)
					else:
						final_r_index = np.append(final_r_index, new_r_index)
						final_r_value = np.append(final_r_value, new_r_value)
					#print(i)
					i = i - 9

	return final_r_value, final_r_index


# opis: sprawdza czasy pomiedzy kolejnymi zalamkami R. Jesli interwal
#       jest mniejszy niz 250 ms zalamek jest uznawany za niepoprawny.
#       Korespondujace z niewlasciwym zalamkiem R granice zespolu QRS
#       sa usuwane z listy.
# nazwa funkcji: FUNC_check_rr_intervals_for_PRINT_ALL
# parametry:
#               r_index - numpy 1D array
#				r_value - numpy 1D array
#               left_T - numpy 1D array
#				right_T - numpy 1D array
#				fs - integer
# zwraca:
#               final_r_value - numpy 1D array
#               final_r_index - numpy 1D array
#               final_left_T - numpy 1D array
#               final_right_T - numpy 1D array
def FUNC_check_rr_intervals_for_PRINT_ALL(r_index, r_value, left_T, right_T, fs):
	r_peaks_to_delete_compute = np.empty(len(r_index))
	del_index = 0
	range_ms = np.ceil(0.36 * fs)

	for i in range(len(r_index)-1):
		rr_interval = r_index[i] - r_index[i - 1]

		if rr_interval < range_ms:
			if r_value[i] < 0.5 * r_value[i - 1]:
				r_peaks_to_delete_compute[del_index] = i
				del_index = del_index + 1

	to_delete = np.empty(del_index)

	for i in range(del_index):
		to_delete[i] = int(r_peaks_to_delete_compute[i])

	final_r_index = np.delete(r_index, to_delete)
	final_r_value = np.delete(r_value, to_delete)

	final_left_T = np.delete(left_T, del_index)
	final_right_T = np.delete(right_T, del_index)

	rr_i = np.empty(len(r_index))

	for i in range(len(r_index)):
		if i > 1:
			rr_i[i] = r_index[i] - r_index[i-1]
			#print(rr_i[i])

			if i > 9:
				avg_rr = (rr_i[i-9]+rr_i[i-8]+rr_i[i-7]+rr_i[i-6]+rr_i[i-5]+rr_i[i-4]+rr_i[i-3]+rr_i[i-2]+rr_i[i-1])/8
				#print(avg_rr)
				if rr_i[i] > 1.66 * avg_rr:
					start = int(r_index[i-1]) + int(0.2 * fs)
					end = int(r_index[i]) - int(0.2 * fs)
					new_r_value = max(ecg_filtered[start:end])
					new_r_index = np.argmax(ecg_filtered[start:end]) + start
					#print(avg_rr)
					#print(end)
					if i < len(final_r_index):
						final_r_index = np.insert(final_r_index, i, new_r_index)
						final_r_value = np.insert(final_r_value, i, new_r_value)
					else:
						final_r_index = np.append(final_r_index, new_r_index)
						final_r_value = np.append(final_r_value, new_r_value)
					#print(i)
					i = i - 9

	return final_r_value, final_r_index, final_left_T, final_right_T


# opis: detekcja zalamka R. Zwraca wartosci i indeksy wyszukanych zalamkow R.
# nazwa funkcji: FUNC_r_detection
# parametry:
#				ecg_filtered - numpy 1D array
#				fs - integer
# zwraca:
#               r_value - numpy 1D array
#               r_index - numpy 1D array
def FUNC_r_detection(ecg_filtered, fs):
	stateDict = {'Done': 1,
				'Signal too short': -2,
				'Incorrect input': -3}

	try:
		if not len(ecg_filtered):
			stateFlag = stateDict['Incorrect input']
			return stateFlag, [[], []]

		izo = FUNC_flat_izo(ecg_filtered)
		ecg_b = FUNC_filter(izo, fs)
		diff = FUNC_diff(ecg_b)
		sqr = FUNC_sqr(diff)
		ecgmf = FUNC_signal_integration(sqr, fs)
		left_T, right_T = FUNC_prepare_for_maxima_search(ecgmf, fs)
		r_value, r_index = FUNC_find_max(left_T, right_T, ecg_filtered)
		final_r_value, final_r_index = FUNC_check_rr_intervals(r_index, r_value, fs, ecg_filtered)

		#print(len(r_index))

		stateFlag = stateDict['Done']

		return final_r_value, final_r_index

	except Exception as e:
		print(f'Module R_PEAKS failed: {e}')
		stateFlag = stateDict['Error']
		return stateFlag, [[], []]


# opis: detekcja zalamka S przy zalozeniu, ze osiaga minimum.
#		Zwraca wartosci i indeksy wyszukanych zalamkow S.
# nazwa funkcji: FUNC_s_detection
# parametry:
#				ecg_filtered - numpy 1D array
#				fs - integer
# zwraca:
#               s_value - numpy 1D array
#               s_index - numpy 1D array
def FUNC_s_detection(ecg_filtered, fs):
	izo = FUNC_flat_izo(ecg_filtered)
	ecg_b = FUNC_filter(izo, fs)
	diff = FUNC_diff(ecg_b)
	sqr = FUNC_sqr(diff)
	ecgmf = FUNC_signal_integration(sqr, fs)
	left_T, right_T = FUNC_prepare_for_maxima_search(ecgmf, fs)
	r_value, r_index = FUNC_find_max(left_T, right_T, ecg_filtered)
	#final_r_value, final_r_index, final_left_T, final_right_T = FUNC_check_rr_intervals_for_PRINT_ALL(r_index, r_value, left_T, right_T, fs)
	s_value, s_index = FUNC_find_min(left_T, right_T, ecg_filtered)

	return s_value, s_index


# opis: Rysuje wykresy poszczegolnych sygnalow z naniesionymi na nie punktami 
# nazwa funkcji: PRINT_all
# parametry:
#				ecg_filtered - numpy 1D array
#				fs - integer
def PRINT_all(ecg_filtered, fs):
	izo = FUNC_flat_izo(ecg_filtered)
	ecg_b = FUNC_filter(izo, fs)
	diff = FUNC_diff(ecg_b)
	sqr = FUNC_sqr(diff)
	ecgmf = FUNC_signal_integration(sqr, fs)
	left_T, right_T = FUNC_prepare_for_maxima_search(ecgmf, fs)
	r_value, r_index = FUNC_find_max(left_T, right_T, ecg_filtered)
	#final_r_value, final_r_index, final_left_T, final_right_T = FUNC_check_rr_intervals_for_PRINT_ALL(r_index, r_value, left_T, right_T, fs)

	print(len(r_index))

	s_value, s_index = FUNC_find_min(left_T, right_T, ecg_filtered)
	qrs_left_values, qrs_right_values, qrs_left_values_ecgmf, qrs_right_values_ecgmf = FUNC_find_qrs_values(left_T, right_T, ecg_filtered, ecgmf)

	plt.figure(1)

	'Zaznaczono przyblizone granice QRS na przefiltrowanym sygnale'
	plt.subplot(611)
	plt.plot(ecg_filtered)
	#plt.plot(left_T,qrs_left_values, marker='o', color='g', ls='')
	#plt.plot(right_T,qrs_right_values, marker='o', color='y', ls='')
	plt.ylabel('ekg_filtered')

	'Wyrownania izolinia'
	plt.subplot(612)
	plt.plot(izo)
	plt.ylabel('izo')

	'Filtracja'
	plt.subplot(613)
	plt.plot(ecg_b)
	plt.ylabel('filt')

	plt.subplot(614)
	plt.plot(diff)
	plt.ylabel('diff')

	plt.subplot(615)
	plt.plot(sqr)
	plt.ylabel('sqr')

	'Zaznaczono przyblizone granice QRS na scalkowanym sygnale'
	plt.subplot(616)
	plt.plot(ecgmf)
	plt.plot(left_T, qrs_left_values_ecgmf, marker='o', color='g', ls='')
	plt.plot(right_T, qrs_right_values_ecgmf, marker='o', color='y', ls='')
	plt.ylabel('ecgmf')

	plt.figure(2)
	plt.plot(ecg_filtered)
	plt.plot(r_index, r_value, marker='x', color='r', ls='')
	plt.plot(s_index, s_value, marker='x', color='b', ls='')
	plt.plot(left_T,qrs_left_values, marker='o', color='g', ls='')
	plt.plot(right_T,qrs_right_values, marker='o', color='y', ls='')
	plt.xlabel('zalamki R, S oraz przyblizone granice QRS')

	plt.figure(3)
	plt.plot(ecg_filtered)
	plt.plot(r_index, r_value, marker='x', color='r', ls='')
	plt.xlabel('zalamki R')

	plt.figure(4)
	plt.plot(FUNC_running_mean(ecgmf, fs))

	plt.show()


# opis: Zaznacza wykryte zalamki R na sygnale EKG
# nazwa funkcji: PRINT_r
# parametry:
#				ecg_filtered - numpy 1D array
#				r_index - numpy 1D array
#				r_value - numpy 1D array
def PRINT_r(ecg_filtered, r_index, r_value):
	plt.figure(5)
	plt.plot(ecg_filtered)
	plt.plot(r_index,r_value, marker='x', color='r', ls='')
	plt.xlabel('zalamki R')

	plt.show()


# opis: Rysuje dowolny przebieg
# nazwa funkcji: PRINT_any
# parametry:
#				signal - numpy 1D array
def PRINT_any(signal):
	plt.figure()
	plt.plot(signal)

	plt.show()


# opis: Rysuje dwa dowolne przebiegi w jednym oknie
# nazwa funkcji: PRINT_any
# parametry:
#				signal - numpy 1D array
#				signal2 - numpy 1D array
def PRINT_compare(signal, signal2):
	plt.figure()
	plt.subplot(211)
	plt.ylabel('ekg')
	plt.plot(signal)
	plt.subplot(212)
	plt.plot(signal2)
	plt.ylabel('filter')

	plt.show()








