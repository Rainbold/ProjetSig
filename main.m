%% Project TS114
%% Computer-aied analys of electrocardiogram signals

%% Vermeulen Gabriel - Peterlin Maxime
%% May 2014

%% 3 ECG visualization

%% 3.1.1 Normal ECG signals
clear all;

ecg_normal1 = load('ecg_normal_1.mat');
ecg_normal2 = load('ecg_normal_2.mat');
ecg_normal3 = load('ecg_normal_3.mat');

subplot 311
N = ecg_normal1.Fs * 4;
plot(ecg_normal1.ecg(1:N))
title('Time evolution of normal ECG signals (1)');
xlabel('Samples (4s)');
ylabel('Amplitude');

subplot 312
N = ecg_normal2.Fs * 4;
plot(ecg_normal2.ecg(1:N))
title('Time evolution of normal ECG signals (2)');
xlabel('Samples (4s)');
ylabel('Amplitude');

subplot 313
N = ecg_normal3.Fs * 4;
plot(ecg_normal3.ecg(1:N))
title('Time evolution of normal ECG signals (3)');
xlabel('Samples (4s)');
ylabel('Amplitude');