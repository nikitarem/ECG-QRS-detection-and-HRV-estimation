clear;
%% начало
disp('========================================================');
disp('Методика цифровой обработки биосигналов сердечного ритма');
disp('Разработал Ремизов Н. В., гр. 6276-120404D');
disp('Проверил Федотов А. А., к.т.н., доцент кафедры ЛБТС');
disp('========================================================');
%загрузка сигнала
signal = load('3.txt')'; %3 вариант

%настройки
Fs = 500; %частота дискретизации
Fmains = 50; %частота сетевого напряжения
lim_plot = [97 100]; %пределы графика детектирования в секундах
win_len = 2; % секунды, ширина окна
min_RR_dist = 0.2; % сек, минимальная дистанция между R-зубцами

%дополнительные расчеты
tm = [0:length(signal) - 1]./Fs; %временнАя ось
Fn = Fs/2; %частота Найквиста
win_len = ceil(win_len * Fs); %отсчетов в окне
min_RR_dist = ceil(min_RR_dist * Fs); %отсчетов между R-зубцами
signal_preproc = signal; %переприсвоение для предобработки

%синтез фильтров
[bm, am] = iirnotch(Fmains/Fn, 0.1); %фильтр сетевой помехи
[bbp, abp] = cheby2(4, 50, [5/Fn 50/Fn], 'bandpass'); %полосовой для выделения QRS
ref_qrs = nikita_create_ref_qrs(Fs); %опорный QRS

%% предобработка и детектирование R-зубцов
% предобработка
signal_preproc = filter(bm, am, signal_preproc); %удаление сетевой помехи
signal_preproc = filtfilt(bbp, abp, signal_preproc); %выделение QRS
signal_preproc = filter(ref_qrs, 1, signal_preproc); %корреляционная фильтрация
signal_preproc = abs(hilbert(signal_preproc)); %преобразование Гильберта

% детектирование R-зубцов
i = 1;
R_locs = [];
while (i < length(signal))
    i = i + win_len;
    
    % формируем скользящее окно
    if ((length(signal) - i) < win_len)
        ecg_detec_qrs = signal_preproc((i - win_len):end);
    else
        ecg_detec_qrs = signal_preproc((i - win_len):i);
    end
    
    % оцениваем порог
    qrs_threshold = 0.4 * max(ecg_detec_qrs);
    
    % ищем qrs
    [~, R_locs_temp] = findpeaks(ecg_detec_qrs, 'MinPeakHeight', qrs_threshold, 'MinPeakDistance', min_RR_dist);
    
    % складируем qrs
    R_locs = [R_locs (R_locs_temp + i - win_len - 1)];
end
disp('Обнаружено R-зубцов');
disp(length(R_locs));

%% анализ ВСР
%подготовка
%расчет NN
NN = diff(R_locs);
%перевод отсчетов в мс
NN = NN * 1/Fs * 1000;
%расчет временнОй оси
tmNN = tm(R_locs);
tmNN = tmNN(2:end);

%удаление аномальных значений
%физиологический критерий
up_thr = 2000; %2000 мс 
down_thr = 200; %200 мс
outlier_counter = sum((NN <= down_thr & NN >= up_thr)); %счетчик выбросов
NN = NN((NN >= down_thr) & (NN <= up_thr));
tmNN = tmNN((NN >= down_thr) & (NN <= up_thr));

%критерий относительного изменения
i = 4;
while (i < length(NN))
   i = i + 1;
   prev4NNmean = mean(NN((i-4):(i-1))); %среднее прошлых 4 отсчетов
   RelDiff = abs(NN(i) - prev4NNmean)/prev4NNmean;
   if (RelDiff > 0.3) %сравниваем с порогом
       NN(i) = [];
       tmNN(i) = [];
       i = i - 1;
       outlier_counter = outlier_counter + 1;
   end
end
disp('Обнаружено выбросов NN');
disp(outlier_counter);

%статистические методы
disp('==============');
disp('Статистисческий анализ ВСР');
SDNN = sqrt((1/length(NN)) * sum((NN - mean(NN)).^2))
RMSSD = sqrt((1/(length(NN) - 1)) * sum(diff(NN).^2))
NN50 = sum(abs(diff(NN)) > 50)
pNN50 = 100 * NN50/length(NN)
CVr = 100 * SDNN/mean(NN)

%геометрические методы
%группируем данные
% NN_binrng - границы столбиков гистограммы
% NN_bincts - число значений в столбиках
NN_binrng = [400:8:1300]; % 8 мс интервал
NN_bincts = histc(NN, NN_binrng); %вычисляем группировку
[max_NN_bincounts, idx_max_NN_bincounts] = max(NN_bincts); %вычисляем индекс максимума

disp('==============');
disp('Геометрический анализ ВСР');
M0 = NN_binrng(idx_max_NN_bincounts) %мода
Am0 = 100 * max_NN_bincounts/(length(NN)) %амплитуда моды
VR = max(NN_binrng(NN_bincts > 0)) - min(NN_binrng(NN_bincts > 0)) %вариационный размах
TrIdx = length(NN)/Am0 %триангулярный индекс ВСР

%спектральные методы
%аппроксимация функции
FsNN = 4; %4 Гц частота дискретизации
i = 1/FsNN; %шаг интерполяции
tmNNaprx = tmNN(1):i:tmNN(end); %временная ось для интерп сигнала
NNaprx = interp1(tmNN, NN, tmNNaprx, 'spline'); %интерполяция

%подготовка к расчету
NNaprx = detrend(NNaprx); %детренд
HamWin = hamming(length(NNaprx)); %окно хэмминга
nfft = 2048; %число точек бпф
df = FsNN/nfft; %шаг по частоте

%расчет спектральной мощности
[psd_estimate, psd_freq] = periodogram(NNaprx, HamWin, nfft, FsNN);

%расчет VLF, LF, HF
Flims = [0.003 0.04 0.15 0.4]; %пределы по частоте
VLF = 0; %исходные значения 
LF = 0;
HF = 0;
i = 0;
f = 0;

while (f <= Flims(4)) %цикл расчета
    f = df * i; %текущая частота
    i = i + 1;
    if ((f >= Flims(1)) && (f < Flims(2)))
       VLF = VLF + psd_estimate(i)*df; %расчет VLF
    elseif ((f >= Flims(2)) && (f < Flims(3)))
       LF = LF + psd_estimate(i)*df; %расчет LF
    elseif (f >= Flims(3))
       HF = HF + psd_estimate(i)*df; %расчет HF
    end 
end

disp('==============');
disp('Непараметрический спектральный анализ ВСР');
VLF = round(VLF)
LF = round(LF)
HF = round(HF)
LFtoHF = LF/HF
TP = HF + LF + VLF

%% графики для отчета
figure('Name','Результаты эксперимента','NumberTitle','off');
uitabgroup(); %группа вкладок
tab1hrv = uitab('Title','Зависимость изменения длительности сердечного цикла от времени');
tab2hrv = uitab('Title','Зависимость изменения СПМ последовательности длительностей сердечного цикла от частоты');
tab1 = uitab('Title','Исходный ЭКС');
tab1_1 = uitab('Title','Опорный сигнал');
tab2 = uitab('Title','Предобработанный ЭКС');
tab3 = uitab('Title','Предобработанный + исходный ЭКС');
tab4 = uitab('Title','ЭКС с детектированными R-зубцами');

% Зависимость изменения длительности сердечного цикла от времени
tab1hrv_axes = axes('parent', tab1hrv);
hold(tab1hrv_axes, 'on'); grid on;
plot(tmNN, NN);
xlim([tmNN(1) 300]); %пределы графика
ylabel('RR-интервалы, мс'); xlabel('Время, с');

% Зависимость изменения СПМ последовательности длительностей сердечного цикла от частоты
tab2hrv_axes = axes('parent', tab2hrv);
hold(tab2hrv_axes, 'on'); grid on;
plot(psd_freq, psd_estimate/10^6);
xlim([0 0.4]); %пределы графика
ylabel('Спектральная плотность мощности, с^2/Гц'); xlabel('Частота, Гц');

%% графики детектора

% исходный ЭКС
tab1_axes = axes('parent', tab1);
hold(tab1_axes, 'on'); grid on;
plot(tm, signal);
xlim(lim_plot); %пределы графика
ylabel('Амплитуда, мВ'); xlabel('Время, с');

% опорный сигнал
tab1_1_axes = axes('parent', tab1_1);
hold(tab1_1_axes, 'on'); grid on;
stem(ref_qrs); plot(ref_qrs, 'color', 'blue');
ylabel('Амплитуда, мВ'); xlabel('Отсчеты');

% предобработанный ЭКС
tab2_axes = axes('parent', tab2);
hold(tab2_axes, 'on'); grid on;
plot(tm, signal_preproc);
xlim(lim_plot);
ylabel('Амплитуда, мВ'); xlabel('Время, с');

% предобработанный + исходный ЭКС
tab3_axes = axes('parent', tab3);
hold(tab3_axes, 'on'); grid on;
plot(tm, signal_preproc);
plot(tm, detrend(signal));
xlim(lim_plot);
ylabel('Амплитуда, мВ'); xlabel('Время, с');

% детектирование
tab4_axes = axes('parent', tab4);
hold(tab4_axes, 'on'); grid on;
plot(tm, signal_preproc);
plot(tm(R_locs), signal_preproc(R_locs), 'ro','MarkerSize',7); %разметка
xlim(lim_plot);
ylabel('Амплитуда, мВ'); xlabel('Время, с');