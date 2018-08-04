import matplotlib.pyplot as plt
import scipy.signal as ss
import numpy as np
import math

window_size = 55
'------------------------------------------------------------------------------'


# opis: oblicza i zwraca pochodna sygnalu.
# nazwa funkcji: FUNC_diff
# parametry:
#               ecg_filtered - numpy 1D array
# zwraca:
#               diff - numpy 1D array
def FUNC_diff(ecg_filtered):
	N = len(ecg_filtered)
	'''POCHODNA---------------------------------------------------------------'''
	#diff = np.diff(ecg_filtered)
	diff = (2*ecg_filtered[5:N]+ecg_filtered[4:N-1]-ecg_filtered[2:N-3]-2*ecg_filtered[1:N-4])/8; 
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
	window = np.ones(window_size)
	'''Calkowanie'''
	temp = ss.lfilter(window,1,sqr)
	ecgmf = ss.medfilt(temp,9)
	ecgmf = ecgmf*dt
	'''Usuniecie opoznienia filtru'''
	delay = math.ceil(len(window) / 2)
	ecgmf = ecgmf[delay:len(ecgmf)]

	return ecgmf



# opis: wyszukuje i zwraca przyblizone granice QRS
#		potrzebne do wyszukania zalamka R (oraz zalamka S).
# nazwa funkcji: FUNC_prepare_for_maxima_search
# parametry:
#               ecgmf - numpy 1D array
# zwraca:
#               left_T - numpy 1D array
#				right_T - numpy 1D array
def FUNC_prepare_for_maxima_search(ecgmf):
	'''Wyszukanie najwyzszej amplitudy'''
	max_A = max(ecgmf)

	'''Budowa tablicy do przeszukiwania'''
	threshold = 0.2
	region_poz = ecgmf>(threshold*max_A)
	region_poz = region_poz*1

	region_poz = np.array([region_poz])

	'''Uzupelnienie zerem'''
	region_poz_LR = np.insert(region_poz, 0, 0)
	region_poz_RL = np.append(region_poz, 0)

	'''SZUKANE MAKSIMOW-------------------------------------------------------'''
	deltaLR = np.diff(region_poz_LR)
	deltaRL = np.diff(region_poz_RL)

	'''Wyznaczenie granic segmentow'''
	left  = np.where(deltaLR==1);
	right = np.where(deltaRL==-1);

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
#               left_T - numpy 1D array
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
		dt = 1/float(fs)
		diff = FUNC_diff(ecg_filtered)
		sqr = FUNC_sqr(diff)
		ecgmf = FUNC_signal_integration(sqr, dt)
		left_T, right_T = FUNC_prepare_for_maxima_search(ecgmf)
		r_value, r_index = FUNC_find_max(left_T, right_T, ecg_filtered)

		stateFlag = stateDict['Done']

		return stateFlag, [r_value, r_index]

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
	dt = 1/float(fs)
	diff = FUNC_diff(ecg_filtered)
	sqr = FUNC_sqr(diff)
	ecgmf = FUNC_signal_integration(sqr, dt)
	left_T, right_T = FUNC_prepare_for_maxima_search(ecgmf)
	s_value, s_index = FUNC_find_min(left_T, right_T, ecg_filtered)

	return s_value, s_index


# opis: Rysuje wykresy poszczegolnych sygnalow z naniesionymi na nie punktami 
# nazwa funkcji: PRINT_all
# parametry:
#				ecg_filtered - numpy 1D array
#				fs - integer
def PRINT_all(ecg_filtered, fs):
	diff = FUNC_diff(ecg_filtered)
	sqr = FUNC_sqr(diff)
	ecgmf = FUNC_signal_integration(sqr, fs)
	left_T, right_T = FUNC_prepare_for_maxima_search(ecgmf)
	r_value, r_index = FUNC_find_max(left_T, right_T, ecg_filtered)
	s_value, s_index = FUNC_find_min(left_T, right_T, ecg_filtered)
	qrs_left_values, qrs_right_values, qrs_left_values_ekgmf, qrs_right_values_ekgmf = FUNC_find_qrs_values(left_T, right_T, ecg_filtered, ecgmf)
	

	plt.figure(1)
	'Zaznaczono przyblizone granice QRS na przefiltrowanym sygnale'
	plt.subplot(411)
	plt.plot(ecg_filtered)
	plt.plot(left_T,qrs_left_values, marker='o', color='g', ls='')
	plt.plot(right_T,qrs_right_values, marker='o', color='y', ls='')
	plt.ylabel('ekg_filtered')

	plt.subplot(412)
	plt.plot(diff)
	plt.ylabel('diff')

	plt.subplot(413)
	plt.plot(sqr)
	plt.ylabel('sqr')

	'Zaznaczono przyblizone granice QRS na scalkowanym sygnale'
	plt.subplot(414)
	plt.plot(ecgmf)
	plt.plot(left_T,qrs_left_values_ekgmf, marker='o', color='g', ls='')
	plt.plot(right_T,qrs_right_values_ekgmf, marker='o', color='y', ls='')
	plt.ylabel('ecgmf')

	plt.figure(2)
	plt.plot(ecg_filtered)
	plt.plot(r_index,r_value, marker='x', color='r', ls='')
	plt.plot(s_index,s_value, marker='x', color='b', ls='')
	plt.plot(left_T,qrs_left_values, marker='o', color='g', ls='')
	plt.plot(right_T,qrs_right_values, marker='o', color='y', ls='')
	plt.xlabel('zalamki R, S oraz przyblizone granice QRS')

	plt.figure(3)
	plt.plot(ecg_filtered)
	plt.plot(r_index,r_value, marker='x', color='r', ls='')
	plt.xlabel('zalamki R')

	plt.show()


# opis: Zaznacza wykryte zalamki R na sygnale EKG
# nazwa funkcji: PRINT_r
# parametry:
#				ecg_filtered - numpy 1D array
#				r_index - numpy 1D array
#				r_value - numpy 1D array
def PRINT_r(ecg_filtered, r_index, r_value):
	plt.figure(4)
	plt.plot(ecg_filtered)
	plt.plot(r_index,r_value, marker='x', color='r', ls='')
	plt.xlabel('zalamki R')

	plt.show()
