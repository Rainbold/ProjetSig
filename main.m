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
title('Time evolution of normal ECG signal (1)');
xlabel('Samples (4s)');
ylabel('Amplitude');
Fs1 = ecg_normal1.Fs

subplot 312
N = ecg_normal2.Fs * 4;
plot(ecg_normal2.ecg(1:N))
title('Time evolution of normal ECG signal (2)');
xlabel('Samples (4s)');
ylabel('Amplitude');
Fs2 = ecg_normal2.Fs

subplot 313
N = ecg_normal3.Fs * 4;
plot(ecg_normal3.ecg(1:N))
title('Time evolution of normal ECG signal (3)');
xlabel('Samples (4s)');
ylabel('Amplitude');
Fs3 = ecg_normal3.Fs

%% 3.1.2 ECG signals with pathologies
ecg_AF = load('ecg_AF.mat');
ecg_VF = load('ecg_VF.mat');
ecg_SSS = load('ecg_SSS.mat');
ecg_PVC = load('ecg_PVC.mat');

subplot 411
N = ecg_AF.Fs * 4;
plot(ecg_AF.ecg(1:N))
title('Time evolution of ECG signal with Atrial Fibrillation');
xlabel('Samples (4s)');
ylabel('Amplitude');
Fs_AF = ecg_AF.Fs

subplot 412
N = ecg_VF.Fs * 4;
plot(ecg_VF.ecg(1:N))
title('Time evolution of ECG signal with Ventricular Fibrillation');
xlabel('Samples (4s)');
ylabel('Amplitude');
Fs_VF = ecg_VF.Fs

subplot 413
N = ecg_SSS.Fs * 4;
plot(ecg_SSS.ecg(1:N))
title('Time evolution of ECG signal with Sick Sinus Syndrome');
xlabel('Samples (4s)');
ylabel('Amplitude');
Fs_SSS = ecg_SSS.Fs

subplot 414
N = ecg_PVC.Fs * 4;
plot(ecg_PVC.ecg(1:N))
title('Time evolution of ECG signal with Premature Ventricular Contraction');
xlabel('Samples (4s)');
ylabel('Amplitude');
Fs_PVC = ecg_PVC.Fs

%% 3.2 Frequency display
%% 3.2.1 Normal ECG signals
S = 15; % Use S seconds of samples

% ECG normal 1
subplot 311
N = Fs1 * S;
dsp_ecg_n1 = abs(fft(ecg_normal1.ecg(1:N))) .^ 2;
plot(([0:N-1]/N-0.5)*Fs1, fftshift(dsp_ecg_n1));
xlim([-Fs1/2,Fs1/2]);

title('Power spectrum of normal ECG signal (1)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% ECG normal 2
subplot 312
N = Fs2 * S; % 4s
dsp_ecg_n2 = abs(fft(ecg_normal2.ecg(1:N))) .^ 2;
plot(([0:N-1]/N-0.5)*Fs2, fftshift(dsp_ecg_n2));
xlim([-Fs2/2,Fs2/2]);

title('Power spectrum of normal ECG signal (2)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% ECG normal 3
subplot 313
N = Fs3 * S; % 4s
dsp_ecg_n3 = abs(fft(ecg_normal3.ecg(1:N))) .^ 2;
plot(([0:N-1]/N-0.5)*Fs3, fftshift(dsp_ecg_n3));
xlim([-Fs3/2,Fs3/2]);

title('Power spectrum of normal ECG signal (3)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

%% 3.2.2 ECG signals with pathologies
S = 50; % Use S seconds of samples

% ECG with AF
subplot 411
N = Fs_AF * S;
dsp_ecg_AF = abs(fft(ecg_AF.ecg(1:N))) .^ 2;
plot(([0:N-1]/N-0.5)*Fs_AF, fftshift(dsp_ecg_AF));
xlim([0;150]);

title('Power spectrum of ECG signal with AF');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% ECG with VF
subplot 412
N = Fs_VF * S; % 4s
dsp_ecg_VF = abs(fft(ecg_VF.ecg(1:N))) .^ 2;
plot(([0:N-1]/N-0.5)*Fs_VF, fftshift(dsp_ecg_VF));
xlim([0;150]);

title('Power spectrum of ECG signal with VF');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot 413
N = Fs_SSS * S; % 4s
dsp_ecg_SSS = abs(fft(ecg_SSS.ecg(1:N))) .^ 2;
plot(([0:N-1]/N-0.5)*Fs_SSS, fftshift(dsp_ecg_SSS));
xlim([0;150]);

title('Power spectrum of ECG signal with SSS');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot 414
N = Fs_PVC * S; % 4s
dsp_ecg_PVC = abs(fft(ecg_PVC.ecg(1:N))) .^ 2;
plot(([0:N-1]/N-0.5)*Fs_PVC, fftshift(dsp_ecg_PVC));
xlim([0;150]);

title('Power spectrum of ECG signal with PVC');
xlabel('Frequency (Hz)');
ylabel('Amplitude');