src/R_PEAKS.py: funkcje szukające załamków R.
FUNC_r_detection(ecg_filtered, fs) zwraca odpowiednio indeksy oraz wartości wykrytych załamków. Jako parametry przyjmuje przefiltrowany sygnał oraz jego częstotliwość (z ECG_BASELINE/ref/out_signalInfo).
Dodatkowo wykrywanie granic QRS oraz załamka S.
Opisy funkcji znajdują się w kodzie.

src/R_PEAKS_test.py: przykład odczytu przefiltrowanych sygnałów, użycia funkcji FUNC_r_detection i rysowania wykresów.