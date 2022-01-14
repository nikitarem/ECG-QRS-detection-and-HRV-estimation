clear;
%% ������
disp('========================================================');
disp('�������� �������� ��������� ����������� ���������� �����');
disp('���������� ������� �. �., ��. 6276-120404D');
disp('�������� ������� �. �., �.�.�., ������ ������� ����');
disp('========================================================');
%�������� �������
signal = load('3.txt')'; %3 �������

%���������
Fs = 500; %������� �������������
Fmains = 50; %������� �������� ����������
lim_plot = [97 100]; %������� ������� �������������� � ��������
win_len = 2; % �������, ������ ����
min_RR_dist = 0.2; % ���, ����������� ��������� ����� R-�������

%�������������� �������
tm = [0:length(signal) - 1]./Fs; %��������� ���
Fn = Fs/2; %������� ���������
win_len = ceil(win_len * Fs); %�������� � ����
min_RR_dist = ceil(min_RR_dist * Fs); %�������� ����� R-�������
signal_preproc = signal; %�������������� ��� �������������

%������ ��������
[bm, am] = iirnotch(Fmains/Fn, 0.1); %������ ������� ������
[bbp, abp] = cheby2(4, 50, [5/Fn 50/Fn], 'bandpass'); %��������� ��� ��������� QRS
ref_qrs = nikita_create_ref_qrs(Fs); %������� QRS

%% ������������� � �������������� R-������
% �������������
signal_preproc = filter(bm, am, signal_preproc); %�������� ������� ������
signal_preproc = filtfilt(bbp, abp, signal_preproc); %��������� QRS
signal_preproc = filter(ref_qrs, 1, signal_preproc); %�������������� ����������
signal_preproc = abs(hilbert(signal_preproc)); %�������������� ���������

% �������������� R-������
i = 1;
R_locs = [];
while (i < length(signal))
    i = i + win_len;
    
    % ��������� ���������� ����
    if ((length(signal) - i) < win_len)
        ecg_detec_qrs = signal_preproc((i - win_len):end);
    else
        ecg_detec_qrs = signal_preproc((i - win_len):i);
    end
    
    % ��������� �����
    qrs_threshold = 0.4 * max(ecg_detec_qrs);
    
    % ���� qrs
    [~, R_locs_temp] = findpeaks(ecg_detec_qrs, 'MinPeakHeight', qrs_threshold, 'MinPeakDistance', min_RR_dist);
    
    % ���������� qrs
    R_locs = [R_locs (R_locs_temp + i - win_len - 1)];
end
disp('���������� R-������');
disp(length(R_locs));

%% ������ ���
%����������
%������ NN
NN = diff(R_locs);
%������� �������� � ��
NN = NN * 1/Fs * 1000;
%������ ��������� ���
tmNN = tm(R_locs);
tmNN = tmNN(2:end);

%�������� ���������� ��������
%��������������� ��������
up_thr = 2000; %2000 �� 
down_thr = 200; %200 ��
outlier_counter = sum((NN <= down_thr & NN >= up_thr)); %������� ��������
NN = NN((NN >= down_thr) & (NN <= up_thr));
tmNN = tmNN((NN >= down_thr) & (NN <= up_thr));

%�������� �������������� ���������
i = 4;
while (i < length(NN))
   i = i + 1;
   prev4NNmean = mean(NN((i-4):(i-1))); %������� ������� 4 ��������
   RelDiff = abs(NN(i) - prev4NNmean)/prev4NNmean;
   if (RelDiff > 0.3) %���������� � �������
       NN(i) = [];
       tmNN(i) = [];
       i = i - 1;
       outlier_counter = outlier_counter + 1;
   end
end
disp('���������� �������� NN');
disp(outlier_counter);

%�������������� ������
disp('==============');
disp('��������������� ������ ���');
SDNN = sqrt((1/length(NN)) * sum((NN - mean(NN)).^2))
RMSSD = sqrt((1/(length(NN) - 1)) * sum(diff(NN).^2))
NN50 = sum(abs(diff(NN)) > 50)
pNN50 = 100 * NN50/length(NN)
CVr = 100 * SDNN/mean(NN)

%�������������� ������
%���������� ������
% NN_binrng - ������� ��������� �����������
% NN_bincts - ����� �������� � ���������
NN_binrng = [400:8:1300]; % 8 �� ��������
NN_bincts = histc(NN, NN_binrng); %��������� �����������
[max_NN_bincounts, idx_max_NN_bincounts] = max(NN_bincts); %��������� ������ ���������

disp('==============');
disp('�������������� ������ ���');
M0 = NN_binrng(idx_max_NN_bincounts) %����
Am0 = 100 * max_NN_bincounts/(length(NN)) %��������� ����
VR = max(NN_binrng(NN_bincts > 0)) - min(NN_binrng(NN_bincts > 0)) %������������ ������
TrIdx = length(NN)/Am0 %������������� ������ ���

%������������ ������
%������������� �������
FsNN = 4; %4 �� ������� �������������
i = 1/FsNN; %��� ������������
tmNNaprx = tmNN(1):i:tmNN(end); %��������� ��� ��� ������ �������
NNaprx = interp1(tmNN, NN, tmNNaprx, 'spline'); %������������

%���������� � �������
NNaprx = detrend(NNaprx); %�������
HamWin = hamming(length(NNaprx)); %���� ��������
nfft = 2048; %����� ����� ���
df = FsNN/nfft; %��� �� �������

%������ ������������ ��������
[psd_estimate, psd_freq] = periodogram(NNaprx, HamWin, nfft, FsNN);

%������ VLF, LF, HF
Flims = [0.003 0.04 0.15 0.4]; %������� �� �������
VLF = 0; %�������� �������� 
LF = 0;
HF = 0;
i = 0;
f = 0;

while (f <= Flims(4)) %���� �������
    f = df * i; %������� �������
    i = i + 1;
    if ((f >= Flims(1)) && (f < Flims(2)))
       VLF = VLF + psd_estimate(i)*df; %������ VLF
    elseif ((f >= Flims(2)) && (f < Flims(3)))
       LF = LF + psd_estimate(i)*df; %������ LF
    elseif (f >= Flims(3))
       HF = HF + psd_estimate(i)*df; %������ HF
    end 
end

disp('==============');
disp('����������������� ������������ ������ ���');
VLF = round(VLF)
LF = round(LF)
HF = round(HF)
LFtoHF = LF/HF
TP = HF + LF + VLF

%% ������� ��� ������
figure('Name','���������� ������������','NumberTitle','off');
uitabgroup(); %������ �������
tab1hrv = uitab('Title','����������� ��������� ������������ ���������� ����� �� �������');
tab2hrv = uitab('Title','����������� ��������� ��� ������������������ ������������� ���������� ����� �� �������');
tab1 = uitab('Title','�������� ���');
tab1_1 = uitab('Title','������� ������');
tab2 = uitab('Title','���������������� ���');
tab3 = uitab('Title','���������������� + �������� ���');
tab4 = uitab('Title','��� � ���������������� R-�������');

% ����������� ��������� ������������ ���������� ����� �� �������
tab1hrv_axes = axes('parent', tab1hrv);
hold(tab1hrv_axes, 'on'); grid on;
plot(tmNN, NN);
xlim([tmNN(1) 300]); %������� �������
ylabel('RR-���������, ��'); xlabel('�����, �');

% ����������� ��������� ��� ������������������ ������������� ���������� ����� �� �������
tab2hrv_axes = axes('parent', tab2hrv);
hold(tab2hrv_axes, 'on'); grid on;
plot(psd_freq, psd_estimate/10^6);
xlim([0 0.4]); %������� �������
ylabel('������������ ��������� ��������, �^2/��'); xlabel('�������, ��');

%% ������� ���������

% �������� ���
tab1_axes = axes('parent', tab1);
hold(tab1_axes, 'on'); grid on;
plot(tm, signal);
xlim(lim_plot); %������� �������
ylabel('���������, ��'); xlabel('�����, �');

% ������� ������
tab1_1_axes = axes('parent', tab1_1);
hold(tab1_1_axes, 'on'); grid on;
stem(ref_qrs); plot(ref_qrs, 'color', 'blue');
ylabel('���������, ��'); xlabel('�������');

% ���������������� ���
tab2_axes = axes('parent', tab2);
hold(tab2_axes, 'on'); grid on;
plot(tm, signal_preproc);
xlim(lim_plot);
ylabel('���������, ��'); xlabel('�����, �');

% ���������������� + �������� ���
tab3_axes = axes('parent', tab3);
hold(tab3_axes, 'on'); grid on;
plot(tm, signal_preproc);
plot(tm, detrend(signal));
xlim(lim_plot);
ylabel('���������, ��'); xlabel('�����, �');

% ��������������
tab4_axes = axes('parent', tab4);
hold(tab4_axes, 'on'); grid on;
plot(tm, signal_preproc);
plot(tm(R_locs), signal_preproc(R_locs), 'ro','MarkerSize',7); %��������
xlim(lim_plot);
ylabel('���������, ��'); xlabel('�����, �');