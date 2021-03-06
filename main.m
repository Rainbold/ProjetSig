%% Project TS114
%% Computer-aied analys of electrocardiogram signals

%% Vermeulen Gabriel - Peterlin Maxime
%% May 2014

%% 3 ECG visualization

%% 3.1.1 Normal ECG signals
clear all; close all;

% Loading ecg normal files
ecg_n1 = load('ecg/ecg_normal_1.mat');
ecg_n2 = load('ecg/ecg_normal_2.mat');
ecg_n3 = load('ecg/ecg_normal_3.mat');

% Sampling frequency of each ecg
Fs1 = ecg_n1.Fs;
Fs2 = ecg_n2.Fs;
Fs3 = ecg_n3.Fs;

% ECG normal 1
subplot 311
N = Fs1 * 4; % number of samples to have 4 seconds
plot((0:N-1)/Fs1, ecg_n1.ecg(1:N));
title('Time evolution of normal ECG signal (1)');
xlabel('Time (s)');
ylabel('Amplitude');

% ECG normal 2
subplot 312
N = Fs2 * 4;
plot((0:N-1)/Fs2, ecg_n2.ecg(1:N))
title('Time evolution of normal ECG signal (2)');
xlabel('Time (s)');
ylabel('Amplitude');

% ECG normal 3
subplot 313
N = Fs3 * 4;
plot((0:N-1)/Fs3, ecg_n3.ecg(1:N))
title('Time evolution of normal ECG signal (3)');
xlabel('Time (s)');
ylabel('Amplitude');

%% 3.1.2 ECG signals with pathologies
close all;

ecg_AF = load('ecg/ecg_AF.mat');
ecg_VF = load('ecg/ecg_VF.mat');
ecg_SSS = load('ecg/ecg_SSS.mat');
ecg_PVC = load('ecg/ecg_PVC.mat');
Fs_AF = ecg_AF.Fs;
Fs_VF = ecg_VF.Fs;
Fs_SSS = ecg_SSS.Fs;
Fs_PVC = ecg_PVC.Fs;

% ECG with AF
subplot 411
N = Fs_AF * 4; % 4 seconds of data
of = Fs_AF * 37; % offset of 37 seconds
plot((of:N-1+of)/Fs_AF, ecg_AF.ecg(1+of:N+of))
title('Time evolution of ECG signal with Atrial Fibrillation');
xlabel('Time (s)');
ylabel('Amplitude');

% ECG with VF
subplot 412
N = Fs_VF * 4;
of = Fs_VF * 280;
plot((of:N-1+of)/Fs_VF, ecg_VF.ecg(1+of:N+of))
title('Time evolution of ECG signal with Ventricular Fibrillation');
xlabel('Time (s)');
ylabel('Amplitude');

% ECG with SSS
subplot 413
N = ecg_SSS.Fs * 4;
of = Fs_SSS * 165;
plot((of:N-1+of)/Fs_SSS, ecg_SSS.ecg(1+of:N+of))
title('Time evolution of ECG signal with Sick Sinus Syndrome');
xlabel('Time (s)');
ylabel('Amplitude');

% ECG with PVC
subplot 414
N = ecg_PVC.Fs * 4;
of = Fs_PVC * 105;
plot((of:N-1+of)/Fs_PVC, ecg_PVC.ecg(1+of:N+of))
title('Time evolution of ECG signal with Premature Ventricular Contraction');
xlabel('Time (s)');
ylabel('Amplitude');

%% 3.2 Frequency display
%% 3.2.1 Normal ECG signals
close all;
s = 15; % Use s seconds of samples

% ECG normal 1
subplot 311
N = Fs1 * s;
dsp_ecg_n1 = abs(fft(ecg_n1.ecg(1+of:N+of))) .^ 2; % Compute the power spectrum
plot(((0:N-1)/N-0.5)*Fs1, fftshift(dsp_ecg_n1));
xlim([-50,50]);

title('Power spectrum of normal ECG signal (1)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% ECG normal 2
subplot 312
N = Fs2 * s;
dsp_ecg_n2 = abs(fft(ecg_n2.ecg(1+of:N+of))) .^ 2;
plot(((0:N-1)/N-0.5)*Fs2, fftshift(dsp_ecg_n2));
xlim([-50,50]);

title('Power spectrum of normal ECG signal (2)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% ECG normal 3
subplot 313
N = Fs3 * s;
dsp_ecg_n3 = abs(fft(ecg_n3.ecg(1+of:N+of))) .^ 2;
plot(((0:N-1)/N-0.5)*Fs3, fftshift(dsp_ecg_n3));
xlim([-50,50]);

title('Power spectrum of normal ECG signal (3)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

%% 3.2.2 ECG signals with pathologies
s = 50; % Use S seconds of samples

% ECG with AF
subplot 221
N = Fs_AF * s;
dsp_ecg_AF = abs(fft(ecg_AF.ecg(1+of:N+of))) .^ 2;
plot(((0:N-1)/N-0.5)*Fs_AF, fftshift(dsp_ecg_AF));
xlim([0;150]);
ylim([0;10^10]);

title('Power spectrum of ECG signal with AF');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% ECG with VF
subplot 222
N = Fs_VF * s;
dsp_ecg_VF = abs(fft(ecg_VF.ecg(1+of:N+of))) .^ 2;
plot(((0:N-1)/N-0.5)*Fs_VF, fftshift(dsp_ecg_VF));
xlim([0;150]);
ylim([0;10^10]);

title('Power spectrum of ECG signal with VF');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% ECG with SSS
subplot 223
N = Fs_SSS * s;
dsp_ecg_SSS = abs(fft(ecg_SSS.ecg(1+of:N+of))) .^ 2;
plot(((0:N-1)/N-0.5)*Fs_SSS, fftshift(dsp_ecg_SSS));
xlim([0;150]);
ylim([0;10^10]);

title('Power spectrum of ECG signal with SSS');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% ECG with PVC
subplot 224
N = Fs_PVC * s; 
dsp_ecg_PVC = abs(fft(ecg_PVC.ecg(1+of:N+of))) .^ 2;
plot(((0:N-1)/N-0.5)*Fs_PVC, fftshift(dsp_ecg_PVC));
xlim([0;150]);
ylim([0;10^10]);

title('Power spectrum of ECG signal with PVC');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

%% 4 Detection of P, QRS and T waves
%% 4.1 QRS
%% 4.1.1 R wave detection
clear all; close all;

ecg_n1 = load('ecg/ecg_normal_1.mat');
ecg_n2 = load('ecg/ecg_normal_2.mat');
ecg_n3 = load('ecg/ecg_normal_3.mat');

Fs1 = ecg_n1.Fs;
Fs2 = ecg_n2.Fs;
Fs3 = ecg_n3.Fs;

%% Method of local maxima

% In order to increase the difference between the R waves and the rest of
% the signal :
ecg1P3 = (ecg_n1.ecg).^3;
ecg2P3 = (ecg_n2.ecg).^3;
ecg3P3 = (ecg_n3.ecg).^3;

s = 8; % number of seconds used
ratio = 30; % percent used to detect peaks

% ECG 1
subplot 311
N = Fs1 * s;
[pks,locs] = FindPeaks(ecg1P3(1:N), ratio, Fs1); % locate the local maxima

plot((0:N-1)/Fs1, ecg1P3(1:N));
hold on
plot((locs-1)/Fs1, pks,'+r');

title('Time evolution of a normal ECG signal by the power of 3 (1)');
legend('ECG signal','R waves');
xlabel('Time (s)'); ylabel('Amplitude');

% ECG 2
subplot 312
N = Fs2 * s;
[pks,locs] = FindPeaks(ecg2P3(1:N), ratio, Fs2); % locate the local maxima

plot((0:N-1)/Fs2, ecg2P3(1:N));
hold on
plot((locs-1)/Fs2, pks,'+r');

title('Time evolution of a normal ECG signal by the power of 3 (2)');
legend('ECG signal','R waves');
xlabel('Time (s)'); ylabel('Amplitude');

% ECG 3
subplot 313
N = Fs3 * s;
[pks,locs] = FindPeaks(ecg3P3(1:N), ratio, Fs3); % locate the local maxima

plot((0:N-1)/Fs3, ecg3P3(1:N));
hold on
plot((locs-1)/Fs3, pks,'+r');

title('Time evolution of a normal ECG signal by the power of 3 (3)');
legend('ECG signal','R waves');
xlabel('Time (s)'); ylabel('Amplitude');

clear ecg1P3;
clear ecg2P3;
clear ecg3P3;
clear locs;
clear pks;
clear N;
clear s;
clear ratio;
%% Method of the derivation
s = 3; % seconds of samples
ratio = 70; % percent used to detect peaks
close all;

figure(1);

% ECG 1
subplot 311
N = Fs1 * s;
[pks,locs] = DerivMeth(ecg_n1.ecg(1:N), ratio, Fs1);

plot((0:N-1)/Fs1, ecg_n1.ecg(1:N));
hold on
plot((locs-1)/Fs1, pks,'+r');

title('Time evolution of a normal ECG signal');
legend('ECG signal','R waves');
xlabel('Time (s)'); ylabel('Amplitude');

% ECG 2
subplot 312
N = Fs2 * s;
[pks,locs] = DerivMeth(ecg_n2.ecg(1:N), ratio, Fs2);

plot((0:N-1)/Fs2, ecg_n2.ecg(1:N));
hold on
plot((locs-1)/Fs2, pks,'+r');

title('Time evolution of a normal ECG signal');
legend('ECG signal','R waves');
xlabel('Time (s)'); ylabel('Amplitude');

% ECG 3
subplot 313
N = Fs3 * s;
[pks,locs] = DerivMeth(ecg_n3.ecg(1:N), ratio, Fs3);

plot((0:N-1)/Fs3, ecg_n3.ecg(1:N));
hold on
plot((locs-1)/Fs3, pks,'+r');

title('Time evolution of a normal ECG signal');
legend('ECG signal','R waves');
xlabel('Time (s)'); ylabel('Amplitude');

clear N;
clear s;
clear locs;
clear pks;
clear ratio;
%% Pan and Tompkins algorithm
%% A)
clear all; close all;
ecg_n1 = load('ecg/ecg_normal_1.mat');
Fs1 = ecg_n1.Fs;

N = 10000;

% Low-pass filter
B_lp = zeros(1,13);
B_lp(1) = 1;
B_lp(7) = -2;
B_lp(13) = 1;

A_lp = [ 1 -2 1 ];

figure(1);
freqz(B_lp, A_lp, N);
title('Low-pass filter amplitude and phase');

% High-pass filter
B_hp = zeros(1,33);
B_hp(1) = -1;
B_hp(17) = 32;
B_hp(18) = -32;
B_hp(33) = 1;

A_hp = [ 1 1 ];

figure(2);
freqz(B_hp, A_hp, N);
title('High-pass filter amplitude and phase');

%% ECG signal filtered
s = 4;
N = Fs1 * s; % s seconds of samples

figure(3);

% ECG normal 1
subplot 231
plot((0:N-1)/N*s, ecg_n1.ecg(1:N));
title('ECG normal signal 1');
xlabel('Time (s)');
ylabel('Amplitude');

% Spectrum of ECG 1
subplot 234
sp_ecg_n1 = abs(fft(ecg_n1.ecg(1:N)));
plot(((0:N-1)/N-0.5)*Fs1, fftshift(sp_ecg_n1));
xlim([-50,50]);
title('Spectrum of ECG normal signal 1');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% Low-pass filter
y = filter(B_lp, A_lp, ecg_n1.ecg(1:N));

subplot 232
plot((0:N-1)/N*s, y);
title('ECG normal signal 1 with low-pass filter');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 235
sp_y = abs(fft(y));
plot(((0:N-1)/N-0.5)*Fs1, fftshift(sp_y));
xlim([-50,50]);
title('Spectrum of ECG with low-pass filter');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% High-pass filter
y1 = filter(B_hp, A_hp, y);

subplot 233
plot((0:N-1)/N*s, y1);
title('ECG normal signal 1 with band-pass filter');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 236
sp_y1 = abs(fft(y1));
plot(((0:N-1)/N-0.5)*Fs1, fftshift(sp_y1));
xlim([-50,50]);
title('Spectrum of ECG with band-pass filter');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

%% B)
clear sp_ecg_n1;
clear A_hp B_hp;
clear A_lp B_lp;

B = [ 1 2 0 -2 -1 ] * Fs1 / 8;
A = 1;

figure(1);
freqz(B,A, 10000);
title('Filter amplitude and phase');

figure(2);

subplot 221
plot((0:N-1)/N*s, y1);
title('ECG normal signal 1 with band-pass filter');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 223
plot(((0:N-1)/N-0.5)*Fs1, fftshift(sp_y1));
xlim([-50,50]);
title('Spectrum of ECG with band-pass filter');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% Filter
y2 = filter(B, A, y1);

subplot 222
plot((0:N-1)/N*s, y2);
title('ECG normal signal 1 filtered');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 224
sp_y2 = abs(fft(y2));
plot(((0:N-1)/N-0.5)*Fs1, fftshift(sp_y2));
xlim([-50,50]);
title('Spectrum of ECG filtered');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

clear A B;
%% C)

y3 = y2 .^ 2;

figure(1);

subplot 221
plot((0:N-1)/N*s, y2);
title('ECG normal signal 1 filtered');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 223
plot(((0:N-1)/N-0.5)*Fs1, fftshift(sp_y2));
xlim([-50,50]);
title('Spectrum of ECG filtered');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot 222
plot((0:N-1)/N*s, y3);
title('ECG normal signal 1 filtered and sqared');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 224
sp_y3 = abs(fft(y3));
plot(((0:N-1)/N-0.5)*Fs1, fftshift(sp_y3));
xlim([-50,50]);
title('Spectrum of ECG filtered');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

%% D)

N2 = 2/10 * Fs1; % Value for the length of the window
Smwi = [];

for n=N2:N;
    t = (1/N2) * sum(y3((n-N2+1):n));
    Smwi = [ Smwi t ];
end

plot((N2:N)/N*s, Smwi);
title('Moving-window integration step');
xlabel('Time (s)');
ylabel('Amplitude');

clear t;
%% E)

max_s = max(Smwi); % Maximum of integrated ECG

threshold = 50; % Threshold used to detect the peaks

[pks, maxima] = FindPeaks(Smwi, threshold, Fs1);

subplot 211
plot((N2:N)/Fs1, Smwi);
hold on;
plot((maxima-1+N2)/Fs1, pks, '+r');
title('Moving-window integration step');
xlabel('Time (s)');
ylabel('Amplitude');


% Detection of the beginin of each peaks
low_threshold = 0.05 * max_s;
data = Smwi;
data(data > low_threshold) = 0;
ds = diff(data);

[pks2, minima] = FindPeaks(ds, 10, Fs1);

minima = minima(1:2:end);
pks2 = pks2(1:2:end);

plot((minima-1+N2)/Fs1, Smwi(minima), '+g');
legend('Integrated ECG','End of rising slopes','Beginning of rising slopes')

% Detection of R peaks
x = [minima maxima]; % groups the minmas/maximas
x = sort(x); % and sort them

R=[];
if(numel(x) ~= 0) % Protection useful for the GUI
    if(x(1) ~= minima(1)) % in case of the first data is not a minima
        I=2:2:(length(x)-1);
    else
       I=1:2:(length(x)-1); % used to get every couple maxima-minima
    end

    delay = N2/2;
    
    for i=I
        y = [ ecg_n1.ecg(x(i)+delay:x(i+1)+delay); x(i)+delay:x(i+1)+delay ]; % store in y the ecg data of the peak and their #
        [maxi, id] = max(y, [], 2); % get the peak and his location with their #
        R = [ R [maxi(1); y(2,id(1))] ]; % store the results
    end
end

subplot 212
plot((0:N-1)/Fs1, ecg_n1.ecg(1:N));
hold on;
plot((R(2,:)-1)/Fs1, R(1,:), '+r');
title('Time evolution of ECG normal 1 signal');
xlabel('Time (s)'); ylabel('Amplitude');
legend('ECG signal', 'R waves');


%% 4.1.2 Q and S waves detection
clear all;
ecg_n1 = load('ecg/ecg_normal_1.mat');
Fs1 = ecg_n1.Fs;

s = 4;
ratio = 50;

N = Fs1 * s;

% R waves detection
[pks,locs] = DerivMeth(ecg_n1.ecg(1:N), ratio, Fs1);

% QS research

d = diff(ecg_n1.ecg(1:N));

M = 25; % Maximum of points to find the 0
Q = [];
R = locs;
S = [];

for i=1:length(R)
    % Q
    for j=1:M
        if(d(R(i) - j) <= 0) % we are looking for a local minima e.g. zero-crossing of the derivate
            Q = [ Q (R(i)-j+1) ];
            break;
        end
    end
    % S
    for j=1:M
        if(d(R(i) + j) >= 0) % same as Q
            S = [ S (R(i)+j) ];
            break;
        end
    end
end

figure(1);
plot((0:N-1)/Fs1, ecg_n1.ecg(1:N));
hold on
plot((Q-1)/Fs1, ecg_n1.ecg(Q),'+g');
plot((R-1)/Fs1, ecg_n1.ecg(R),'+r');
plot((S-1)/Fs1, ecg_n1.ecg(S),'+y');
text(Q(1)/Fs1,ecg_n1.ecg(Q(1))*1.1,'Q'), text(R(1)/Fs1,ecg_n1.ecg(R(1))*1.03,'R'), text(S(1)/Fs1,ecg_n1.ecg(S(1))*1.1,'S');
title('Time evolution of a normal ECG signal');
legend('ECG signal','Q waves', 'R waves', 'S waves');
xlabel('Time (s)'); ylabel('Amplitude');

%% QRS research for the 3 normal ECG signals
clear all;
ecg_n1 = load('ecg/ecg_normal_1.mat');
ecg_n2 = load('ecg/ecg_normal_2.mat');
ecg_n3 = load('ecg/ecg_normal_3.mat');

Fs1 = ecg_n1.Fs;
Fs2 = ecg_n2.Fs;
Fs3 = ecg_n3.Fs;

figure(2);

s = 10; % s seconds of samples
ratio = 30; % percent used to detect peak
ratio_deviation = 1/12; % Fs * ratio_deviation : maximum deviation in point, to detect Q and S from R

% ECG 1
subplot 311
N = Fs1 * s;

[Q, R, S] = QRS(2, ecg_n1.ecg(1:N), Fs1, ratio, ratio_deviation);

plot((0:N-1)/Fs1, ecg_n1.ecg(1:N));
hold on;
plot((Q-1)/Fs1, ecg_n1.ecg(Q),'+g');
plot((R-1)/Fs1, ecg_n1.ecg(R),'+r');
plot((S-1)/Fs1, ecg_n1.ecg(S),'+k');
title('Time evolution of normal ECG signal (1)');
xlabel('Time (s)');
ylabel('Amplitude');
text(Q(1)/Fs1,ecg_n1.ecg(Q(1))*1.1,'Q'), text(R(1)/Fs1,ecg_n1.ecg(R(1))*1.03,'R'), text(S(1)/Fs1,ecg_n1.ecg(S(1))*1.1,'S')

% ECG 2
subplot 312
N = Fs2 * s;

[Q, R, S] = QRS(2, ecg_n2.ecg(1:N), Fs2, ratio, ratio_deviation);

plot((0:N-1)/Fs2, ecg_n2.ecg(1:N))
hold on;
plot((Q-1)/Fs2, ecg_n2.ecg(Q),'+g');
plot((R-1)/Fs2, ecg_n2.ecg(R),'+r');
plot((S-1)/Fs2, ecg_n2.ecg(S),'+k');
title('Time evolution of normal ECG signal (2)');
xlabel('Time (s)');
ylabel('Amplitude');
text(Q(1)/Fs2,ecg_n2.ecg(Q(1))*1.1,'Q'), text(R(1)/Fs2,ecg_n2.ecg(R(1))*1.03,'R'), text(S(1)/Fs2,ecg_n2.ecg(S(1))*1.1,'S')

% ECG 3
subplot 313
N = Fs3 * s;

[Q, R, S] = QRS(2, ecg_n3.ecg(1:N), Fs3, ratio, ratio_deviation);

plot((0:N-1)/Fs3, ecg_n3.ecg(1:N))
hold on;
plot((Q-1)/Fs3, ecg_n3.ecg(Q),'+g');
plot((R-1)/Fs3, ecg_n3.ecg(R),'+r');
plot((S-1)/Fs3, ecg_n3.ecg(S),'+k');
title('Time evolution of normal ECG signal (3)');
xlabel('Time (s)');
ylabel('Amplitude');
text(Q(1)/Fs3,ecg_n3.ecg(Q(1))*1.1,'Q'), text(R(1)/Fs3,ecg_n3.ecg(R(1))*1.03,'R'), text(S(1)/Fs3,ecg_n3.ecg(S(1))*1.1,'S')

%% 4.2 P and T waves detection
%% 4.2.1 Filters
clear all;
ecg_n1 = load('ecg/ecg_normal_1.mat');
Fs1 = ecg_n1.Fs;

% Differiator
B_g1 = [ 1 0 0 0 0 0 -1];
A_g1 = 1;

figure(1);
freqz(B_g1, A_g1);
title('Amplitude ans phase of the first filter');

% Low-pass filter
B_g2 = [ 1 0 0 0 0 0 0 0 -1 ];
A_g2 = [ 1 -1 ];

figure(2);
freqz(B_g2, A_g2);
title('Amplitude ans phase of the second filter');

%% 4.2.2 Filtered data
S = 4;
N = S * Fs1;
ratio = 40;

[~,R] = DerivMeth(ecg_n1.ecg(1:N), ratio, Fs1);

% We take 2 R as an interval
R1 = R(1);
R2 = R(2);

E = R2 - R1;
N2 = floor(E * 0.7);


N1 = R1+(Fs1/20);
data = ecg_n1.ecg(N1:N2+R1);

y = filter(B_g2, A_g2, data);
y = y(7:end);
y1 = filter(B_g1, A_g1, y);
y1 = y1(9:end);

delay = (9+7)/2;

figure(3);
subplot 211
plot( (N1+delay:(length(y1)+N1+delay-1))/Fs1, y1);
title('Time evolution of normal ECG signal 1 filtered');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 212
N = 1000;
sp = abs(fft(y1, N));
plot(((0:N-1)/N-0.5)*Fs1, fftshift(sp));
%xlim([-50,50]);
title('Spectrum of ECG normal signal 1 filtered');
xlabel('Frequency (Hz)');
ylabel('Amplitude')

%% Detection of the T waves
figure(1);

subplot 212
plot((N1+delay:N1+delay+length(y1)-1)/Fs1,y1);
hold on;
title('Time evolution of a normal ECG signal');
xlabel('Time (s)'); ylabel('Amplitude');

d = [y1];
[M, idM] = max(d);
[m, idm] = min(d);
plot((N1+delay+idM-1)/Fs1, M, 'r+');
plot((N1+delay+idm-1)/Fs1, m, 'r+');

p = y1(idM:idm);
plot((N1+delay+idM-1:N1+delay+idm-1)/Fs1,p, 'r');
title('Time evolution of a normal ECG signal filtered (1)');
xlabel('Time (s)'); ylabel('Amplitude');

%find zero crossings
t1=p(1:end-1);
t2=p(2:end);
tt=t1.*t2;
indx=find(tt<0);
indx = indx(1);

plot((N1+delay+indx+idM-2)/Fs1, p(indx),'+g');
N = length(data);

gd = 6;

subplot 211
plot((N1+delay:N1+delay+N-1)/Fs1,data);
hold on
plot((N1+delay+indx+idM +gd -1)/Fs1, data(indx+idM+gd), '+g');
title('Time evolution of a normal ECG signal');
legend('ECG signal','T waves');
xlabel('Time (s)'); ylabel('Amplitude');
%% P wave
figure(1);

data = ecg_n1.ecg((R1+N2):(R2-(Fs1/20)));

y = filter(B_g2, A_g2, data);
y = y(7:end);
y1 = filter(B_g1, A_g1, y);
y1 = y1(9:end);

delay = (9+7)/2;

subplot 313
plot((((R1+N2+7):R2-7-(Fs1/20))-1)/Fs1,y1);
hold on
title('Time evolution of a normal ECG signal');
xlabel('Time (s)'); ylabel('Amplitude');

d = [y1; ((R1+N2+7):R2-7-(Fs1/20))];
[M, idM] = max(d, [], 2);
[m, idm] = min(d, [], 2);
plot((d(2,idM(1))-1)/Fs1, M(1), 'r+');
plot((d(2,idm(1))-1)/Fs1, m(1), 'r+');

p = y1(idM(1):idm(1));

plot(((d(2,idM(1)):d(2,idm(1)))-1)/Fs1,p, 'r');
title('Time evolution of a normal ECG signal filtered #1 (y1)');
xlabel('Time (s)'); ylabel('Amplitude');

%find zero crossings
t1=p(1:end-1);
t2=p(2:end);
tt=t1.*t2;
indx=find(tt<0);
indx = indx(1);

plot(((d(2,idM(1)))+indx-2)/Fs1, p(indx),'+g');
N = length(data);

gd = 6;

subplot 311
plot(((R1+N2):(R1+N2+N)-1)/Fs1,data);
hold on
plot((R1+N2+indx+gd-1)/Fs1, data(indx+gd), '+g');
title('Time evolution of a normal ECG signal');
legend('ECG signal','T waves');
xlabel('Time (s)'); ylabel('Amplitude');

subplot 312
plot(((R1+N2):(R1+N2+length(y))-1)/Fs1,y);
hold on
plot((R1+N2+indx+gd-1)/Fs1,y(indx+gd), '+g');
title('Time evolution of the EG signal #1 filtered (y)');
legend('ECG signal','T waves');
xlabel('Time (s)'); ylabel('Amplitude');

%% PQRST display of the 3 normal ECG signals
clear all;

ecg_n1 = load('ecg/ecg_normal_1.mat');
ecg_n2 = load('ecg/ecg_normal_2.mat');
ecg_n3 = load('ecg/ecg_normal_3.mat');

Fs1 = ecg_n1.Fs;
Fs2 = ecg_n2.Fs;
Fs3 = ecg_n3.Fs;

s = 4; % S seconds of samples
ratio = 40;
ratio_d = 1/12; % Fs * ratio_deviation : maximum deviation in point, to detect Q and S from R


% ECG 1
subplot 311
N = Fs1 * s;

[ P, Q, R, S, T ] = PQRST(2, ecg_n1.ecg(1:N), Fs1, ratio, ratio_d);

plot((0:N-1)/Fs1, ecg_n1.ecg(1:N));
hold on;
plot((P-1)/Fs1, ecg_n1.ecg(P),'+y');
plot((Q-1)/Fs1, ecg_n1.ecg(Q),'+g');
plot((R-1)/Fs1, ecg_n1.ecg(R),'+r');
plot((S-1)/Fs1, ecg_n1.ecg(S),'+k');
plot((T-1)/Fs1, ecg_n1.ecg(T),'+r');

title('Time evolution of normal ECG signal (1)');
xlabel('Time (s)');
ylabel('Amplitude');
text(P(1)/Fs1,ecg_n1.ecg(P(1))*5,'P');
text(Q(1)/Fs1,ecg_n1.ecg(Q(1))*1.5,'Q');
text(R(1)/Fs1,ecg_n1.ecg(R(1))*1.1,'R');
text(S(1)/Fs1,ecg_n1.ecg(S(1))*1.4,'S');
text(T(1)/Fs1,ecg_n1.ecg(T(1))*1.5,'T');

% ECG 1
subplot 312
N = Fs2 * s;

[ P, Q, R, S, T ] = PQRST(2, ecg_n2.ecg(1:N), Fs2, ratio, ratio_d);

plot((0:N-1)/Fs2, ecg_n2.ecg(1:N));
hold on;
plot((P-1)/Fs2, ecg_n2.ecg(P),'+y');
plot((Q-1)/Fs2, ecg_n2.ecg(Q),'+g');
plot((R-1)/Fs2, ecg_n2.ecg(R),'+r');
plot((S-1)/Fs2, ecg_n2.ecg(S),'+k');
plot((T-1)/Fs2, ecg_n2.ecg(T),'+r');

title('Time evolution of normal ECG signal (2)');
xlabel('Time (s)');
ylabel('Amplitude');
text(P(1)/Fs2,ecg_n2.ecg(P(1))*5,'P');
text(Q(1)/Fs2,ecg_n2.ecg(Q(1))*1.5,'Q');
text(R(1)/Fs2,ecg_n2.ecg(R(1))*1.1,'R');
text(S(1)/Fs2,ecg_n2.ecg(S(1))*1.4,'S');
text(T(1)/Fs2,ecg_n2.ecg(T(1))*1.5,'T');

% ECG 3
subplot 313
N = Fs3 * s;
offset = Fs3 * 2;

[ P, Q, R, S, T ] = PQRST(2, ecg_n3.ecg(1+offset:N+offset), Fs3, ratio, ratio_d);

plot((0+offset:N-1+offset)/Fs3, ecg_n3.ecg(1+offset:N+offset));
hold on;
plot((P-1+offset)/Fs3, ecg_n3.ecg(P+offset),'+y');
plot((Q-1+offset)/Fs3, ecg_n3.ecg(Q+offset),'+g');
plot((R-1+offset)/Fs3, ecg_n3.ecg(R+offset),'+r');
plot((S-1+offset)/Fs3, ecg_n3.ecg(S+offset),'+k');
plot((T-1+offset)/Fs3, ecg_n3.ecg(T+offset),'+r');

title('Time evolution of normal ECG signal (3)');
xlabel('Time (s)');
ylabel('Amplitude');
text((P(1)+offset)/Fs3,ecg_n3.ecg(P(2)+offset)*-10,'P');
text((Q(1)+offset-40)/Fs3,ecg_n3.ecg(Q(2)+offset)*1.5,'Q');
text((R(1)+offset)/Fs3,ecg_n3.ecg(R(2)+offset)*1.1,'R');
text((S(1)+offset+10)/Fs3,ecg_n3.ecg(S(2)+offset)*0.8,'S');
%text((T(1)+offset)/Fs3,ecg_n3.ecg(T(2)+offset)*1.5,'T');
% It can't find T

%% 5 Automatic identification of cardiac pathologies
%% 5.1 Tachycardia / Bradycardia
clear all;

ecg_n1 = load('ecg/ecg_normal_1.mat');
ecg_n2 = load('ecg/ecg_normal_2.mat');
ecg_n3 = load('ecg/ecg_normal_3.mat');
ecg_AF = load('ecg/ecg_AF.mat');
ecg_VF = load('ecg/ecg_VF.mat');
ecg_SSS = load('ecg/ecg_SSS.mat');
ecg_PVC = load('ecg/ecg_PVC.mat');

Fs1 = ecg_n1.Fs;
Fs2 = ecg_n2.Fs;
Fs3 = ecg_n3.Fs;
Fs_AF = ecg_AF.Fs;
Fs_VF = ecg_VF.Fs;
Fs_SSS = ecg_SSS.Fs;
Fs_PVC = ecg_PVC.Fs;

ratio = 40;
ratio_d = 1/12; % Fs * ratio_deviation : maximum deviation in point, to detect Q and S from R
s = 15;

% ECG normal 1
figure(1);
N = s * Fs1;
[ ~, R, ~] = QRS(2, ecg_n1.ecg(1:N), Fs1, ratio, ratio_d);
N = length(R) - 1;
Ds = diff(R);
Dbpm = (60 * Fs1) * Ds.^-1;
subplot 211
plot((R(1:N)/R(N)*s), Ds/Fs1, '+-');
xlabel('Time (s)'); ylabel('seconds');
title('Cardiac rhythm of ECG normal 1');
subplot 212
plot((R(1:N)/R(N)*s), Dbpm, '+r-');
xlabel('Time (s)'); ylabel('beats/minute');
title('Cardiac rhythm of ECG normal 1');
ylim([0 100]);

% ECG normal 2
figure(2);
N = s * Fs2;
[ ~, R, ~] = QRS(2, ecg_n2.ecg(1:N), Fs2, ratio, ratio_d);
N = length(R) - 1;
Ds = diff(R);
Dbpm = (60 * Fs2) * Ds.^-1;
subplot 211
plot(R(1:N)/R(N)*s, Ds/Fs2, '+-');
xlabel('Time (s)'); ylabel('seconds');
title('Cardiac rhythm of ECG normal 2');
subplot 212
plot((R(1:N)/R(N)*s), Dbpm, '+r-');
xlabel('Time (s)'); ylabel('beats/minute');
title('Cardiac rhythm of ECG normal 2');
ylim([0 100]);

% ECG with VF
figure(3)
N = s * Fs_VF;
[ ~, R, ~] = QRS(2, ecg_VF.ecg(1:N), Fs_VF, 45, ratio_d);
N = length(R) - 1;
Ds = diff(R);
Dbpm = (60 * Fs_VF) * Ds.^-1;
subplot 211
plot((R(1:N)/R(N)*s), Ds/Fs_VF, '+-');
xlabel('Time (s)'); ylabel('seconds');
title('Cardiac rhythm of ECG with VF');
subplot 212
plot((R(1:N)/R(N)*s), Dbpm, '+r-');
xlabel('Time (s)'); ylabel('beats/minute');
title('Cardac rhythm of ECG with VF');
ylim([0 100]);

% ECG with SSS
figure(4)
N = s * Fs_SSS;
[ ~, R, ~] = QRS(2, ecg_SSS.ecg(1:N), Fs_SSS, 45, ratio_d);
N = length(R) - 1;
Ds = diff(R);
Dbpm = (60 * Fs_SSS) * Ds.^-1;
subplot 211
plot((R(1:N)/R(N)*s), Ds/Fs_VF, '+-');
xlabel('Time (s)'); ylabel('seconds');
title('Cardiac rhythm of ECG with SSS');
subplot 212
plot((R(1:N)/R(N)*s), Dbpm, '+r-');
xlabel('Time (s)'); ylabel('beats/minute');
title('Cardiac rhythm of ECG with SSS');
ylim([0 100]);

clear A B N Q R S ratio ratio_d s Ds Dbpm;
%% 5.2 Heart rate variability
close all;
s = 22;
ratio = 40;
ratio_d = 1/12;
Fs = Fs1;

N = s * Fs;

[StdDev, dr, v, bpm] = HRV(1, ecg_n2.ecg(1:N), Fs, ratio, ratio_d);

% Standar Deviation
StdDev

% Histogram of occurences
figure(1)
hist(dr);
title('Histogram of occurences');

figure(1);
subplot 211
plot(v);
title('Process v(t)');

% Spectrum
N = 10^5;
subplot 212;
sp = fftshift(abs(fft((v-mean(v)),N)));
plot(((0:N-1)/N-0.5)*Fs, sp);
xlim([0 0.4]);
title('Spectrum of v(t)');
xlabel('Frequency (Hz)'); ylabel('Amplitude');
hold on;
plot(bpm/60, sp(N/2 +1+ (bpm*N)/(60*Fs)), '+r');

% Respiratory rhythm
bpm

%bpm = Respiratory_rhythm(R, Fs)

%% 5.3 Ectopic beat
s = 18;
ratio = 40;

N = Fs_PVC * s;
[pks,R] = DerivMeth(ecg_PVC.ecg(1:N), ratio, Fs_PVC);

figure(1)
subplot 211
plot((0:N-1)/N*s, ecg_PVC.ecg(1:N));
title('Time evolution of ECG with PVC');
xlabel('Time (s)'); ylabel('Amplitude');
hold on

% Le threshold ideal serait de l'ordre d'une variation max de 5 ou 10 bpm
% une variation de 5 bpm <=> 13 ou 14 point
% Mettons alors 15 points pour essayer

dr = diff(R);
ddr = abs(diff(dr));
subplot 212
plot(R(3:end)/R(end)*s,ddr)
hold on
threshold = 15;
ectopic = find(ddr > threshold);
delay = 2;
plot(R(ectopic+delay)/R(end)*s, ddr(ectopic), '+r');
title('Second derivative of R waves');
xlabel(''); ylabel('');


subplot 211
plot((R(ectopic(1)+delay):R(ectopic(1)+delay+1))/N*s, ecg_PVC.ecg(R(ectopic(1)+delay):R(ectopic(1)+delay+1)), 'r');
plot((R(ectopic(2)+delay):R(ectopic(2)+delay+1))/N*s, ecg_PVC.ecg(R(ectopic(2)+delay):R(ectopic(2)+delay+1)), 'r');


%% 5.4 Fibrillation
sec = 10;
close all;

% Autocovariance : ECG normal 1
N = sec * Fs1;
of = 4 * Fs1; % Offset
ecg_data = ecg_n1.ecg(1+of:N+of);
[~, R, ~] = QRS(1, ecg_data, Fs1, 40, 1/12);
gamma = [];
dr = diff(R);
drm = mean(dr);
s = 0;
N2 = length(dr);
for k=1:N2-2
    for n=1:(N2-k-1)
       s = s + ( (dr(n+k) - drm).*(dr(n) - drm) );
    end
   gamma = [ gamma (1/(N2-k-1) * s) ];
end

subplot 211
plot(R(1:length(gamma)), gamma, 'r');
title('Autocovariance of ECG normal 1 signal');
subplot 212
plot(ecg_data);
title('ECG signal #1');
hold on
plot(R,ecg_data(R), '+r');
xlabel('Samples'); ylabel('Amplitude');
legend('ECG', 'R waves');

% Autocovariance : ECG with AF
figure
N = sec * Fs_AF;
of = 50 * Fs_AF; % Offset
ecg_data = ecg_AF.ecg(1+of:N+of);
[~, R, ~] = QRS(3, ecg_data, Fs_AF, 30, 1/12);
gamma = [];
dr = diff(R);
drm = mean(dr);
s = 0;
N2 = length(dr);
for k=1:N2-2
    for n=1:(N2-k-1)
       s = s + ( (dr(n+k) - drm).*(dr(n) - drm) );
    end
   gamma = [ gamma (1/(N2-k-1) * s) ];
end

subplot 211
plot(R(1:length(gamma)), gamma, 'r');
title('Autocovariance of ECG normal 1 signal');

subplot 212
plot(ecg_data);
hold on
plot(R,ecg_data(R), '+r');
xlabel('Samples'); ylabel('Amplitude');
legend('ECG', 'R waves');


%% No P with AF
s = 4;
N = Fs_AF * s;
ratio = 40;
ratio_d = 1/12;

[ P, Q, R, S, T ] = PQRST(2, ecg_AF.ecg(1:N), Fs_AF, ratio, ratio_d);

plot((0:N-1)/Fs_AF, ecg_AF.ecg(1:N));
hold on;
plot((P-1)/Fs_AF, ecg_AF.ecg(P),'+y');
plot((Q-1)/Fs_AF, ecg_AF.ecg(Q),'+g');
plot((R-1)/Fs_AF, ecg_AF.ecg(R),'+r');
plot((S-1)/Fs_AF, ecg_AF.ecg(S),'+k');
plot((T-1)/Fs_AF, ecg_AF.ecg(T),'+r');

title('Time evolution of normal ECG signal (1)');
xlabel('Time (s)');
ylabel('Amplitude');
%text(P(1)/Fs_AF, ecg_AF.ecg(P(1))*5,'P');
text(Q(1)/Fs_AF-0.03, ecg_AF.ecg(Q(1))*1.1,'Q');
text(R(1)/Fs_AF+0.05, ecg_AF.ecg(R(1)),'R');
text(S(1)/Fs_AF, ecg_AF.ecg(S(1))*1.1,'S');
text(T(1)/Fs_AF, ecg_AF.ecg(T(1))*1.5,'T');

P

%% Ventricular fibrilation
clear all; close all;

ecg_VF = load('ecg/ecg_VF.mat');
Fs = ecg_VF.Fs;

D = 225 * Fs;
F = 235 * Fs;
ecg = ecg_VF.ecg(D:F);

subplot 211
plot(((D:F)-1)/Fs, ecg);
title('Time evolution of an ECG with VF');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 212
N = length(ecg);
sp_ecg = abs(fft(ecg));
sp_ecg(1) = 0;
sp_ecg = fftshift(sp_ecg);
plot(((0:N-1)/N-0.5)*Fs, sp_ecg);

N = length(sp_ecg);
sp = sp_ecg(floor(N/2):end);
id_max_sp = find(sp == max(sp))
title('Power spectrum of normal ECG signal (2)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

rhythm = id_max_sp * Fs / N * 60

%% 6 ECG denoising
clear all;
ecg_PL_file = load('ecg/ecg_noisePL.mat');

Fs = ecg_PL_file.Fs;

s = 4;
N = s * Fs;

ecg = ecg_PL_file.ecg(1:N);


subplot 211
plot((0:N-1)/Fs, ecg);
title('Time evolution of ECG signal');
xlabel('Time (s)');
ylabel('Amplitude');

N = 1000000;
subplot 212
plot(((0:N-1)/N-0.5)*Fs, fftshift(abs(fft(ecg, N)).^2));
ylim([0 1000000]);
title('Spectrum of ECG signal');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% sinus f = 50Hz

%% Filter
[b, a] = butter(2, 0.02, 'low');
y = filter(b, a, ecg);
N = s * Fs;

subplot 221
plot((0:N-1)/Fs, ecg);
title('Time evolution of ECG signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 222
plot((0:N-1)/Fs, y);
title('Time evolution of ECG sgnal filtered');
xlabel('Time (s)');
ylabel('Amplitude');

N = 10^6;

subplot 223
plot(((0:N-1)/N-0.5)*Fs, fftshift(abs(fft(ecg, N)).^2));
ylim([0 1000000]);
title('Spectrum of ECG signal');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot 224
plot(((0:N-1)/N-0.5)*Fs, fftshift(abs(fft(y, N)).^2));
ylim([0 1000000]);
title('Spectrum of ECG signal filtered');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

%% noiseBL
clear all;
ecg_BL_file = load('ecg/ecg_noiseBL.mat');

Fs = ecg_BL_file.Fs;

s = 4;
N = s * Fs;

ecg = ecg_BL_file.ecg(1:N);


subplot 211
plot((0:N-1)/Fs, ecg);
title('Time evolution of ECG signal');
xlabel('Time (s)');
ylabel('Amplitude');

N = 1000000;
subplot 212
plot(((0:N-1)/N-0.5)*Fs, fftshift(abs(fft(ecg, N)).^2));
ylim([0 1000000]);
title('Spectrum of ECG signal');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

%% Filter
[b, a] = cheby1(2, 0.1, 0.01, 'high');

y = filter(b, a, ecg);

N = s * Fs;

subplot 221
plot((0:N-1)/Fs, ecg);
title('Time evolution of ECG signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 222
plot((0:N-1)/Fs, y);
title('Time evolution of ECG sgnal filtered');
xlabel('Time (s)');
ylabel('Amplitude');

N = 10^6;

subplot 223
plot(((0:N-1)/N-0.5)*Fs, fftshift(abs(fft(ecg, N)).^2));
ylim([0 1000000]);
title('Spectrum of ECG signal');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot 224
plot(((0:N-1)/N-0.5)*Fs, fftshift(abs(fft(y, N)).^2));
ylim([0 1000000]);
title('Spectrum of ECG signal filtered');
xlabel('Frequency (Hz)');
ylabel('Amplitude');


%% 7 Graphical interface

% See GUI.m and GUI.fig
% enter GUI to start the GUI



