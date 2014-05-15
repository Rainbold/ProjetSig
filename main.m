%% Project TS114
%% Computer-aied analys of electrocardiogram signals

%% Vermeulen Gabriel - Peterlin Maxime
%% May 2014

%% 3 ECG visualization

%% 3.1.1 Normal ECG signals
clear all;

ecg_n1 = load('ecg/ecg_normal_1.mat');
ecg_n2 = load('ecg/ecg_normal_2.mat');
ecg_n3 = load('ecg/ecg_normal_3.mat');

Fs1 = ecg_n1.Fs;
Fs2 = ecg_n2.Fs;
Fs3 = ecg_n3.Fs;


subplot 311
N = Fs1 * 4;
plot((0:N-1)/Fs1, ecg_n1.ecg(1:N));
title('Time evolution of normal ECG signal (1)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 312
N = Fs2 * 4;
plot((0:N-1)/Fs2, ecg_n2.ecg(1:N))
title('Time evolution of normal ECG signal (2)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 313
N = Fs3 * 4;
plot((0:N-1)/Fs3, ecg_n3.ecg(1:N))
title('Time evolution of normal ECG signal (3)');
xlabel('Time (s)');
ylabel('Amplitude');

%% 3.1.2 ECG signals with pathologies
ecg_AF = load('ecg/ecg_AF.mat');
ecg_VF = load('ecg/ecg_VF.mat');
ecg_SSS = load('ecg/ecg_SSS.mat');
ecg_PVC = load('ecg/ecg_PVC.mat');
Fs_AF = ecg_AF.Fs;
Fs_VF = ecg_VF.Fs;
Fs_SSS = ecg_SSS.Fs;
Fs_PVC = ecg_PVC.Fs;

subplot 411
N = Fs_AF * 4;
plot((0:N-1)/Fs_AF, ecg_AF.ecg(1:N))
title('Time evolution of ECG signal with Atrial Fibrillation');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 412
N = Fs_VF * 4;
plot((0:N-1)/Fs_VF, ecg_VF.ecg(1:N))
title('Time evolution of ECG signal with Ventricular Fibrillation');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 413
N = ecg_SSS.Fs * 4;
plot((0:N-1)/Fs_SSS, ecg_SSS.ecg(1:N))
title('Time evolution of ECG signal with Sick Sinus Syndrome');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 414
N = ecg_PVC.Fs * 4;
plot((0:N-1)/Fs_PVC, ecg_PVC.ecg(1:N))
title('Time evolution of ECG signal with Premature Ventricular Contraction');
xlabel('Time (s)');
ylabel('Amplitude');

%% 3.2 Frequency display
%% 3.2.1 Normal ECG signals
s = 15; % Use S seconds of samples

% ECG normal 1
subplot 311
N = Fs1 * s;
dsp_ecg_n1 = abs(fft(ecg_n1.ecg(1:N))) .^ 2;
plot(((0:N-1)/N-0.5)*Fs1, fftshift(dsp_ecg_n1));
xlim([-50,50]);

title('Power spectrum of normal ECG signal (1)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% ECG normal 2
subplot 312
N = Fs2 * s; % 4s
dsp_ecg_n2 = abs(fft(ecg_n2.ecg(1:N))) .^ 2;
plot(((0:N-1)/N-0.5)*Fs2, fftshift(dsp_ecg_n2));
xlim([-50,50]);

title('Power spectrum of normal ECG signal (2)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% ECG normal 3
subplot 313
N = Fs3 * s; % 4s
dsp_ecg_n3 = abs(fft(ecg_n3.ecg(1:N))) .^ 2;
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
dsp_ecg_AF = abs(fft(ecg_AF.ecg(1:N))) .^ 2;
plot(((0:N-1)/N-0.5)*Fs_AF, fftshift(dsp_ecg_AF));
xlim([0;150]);
ylim([0;10^10]);

title('Power spectrum of ECG signal with AF');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% ECG with VF
subplot 222
N = Fs_VF * s; % 4s
dsp_ecg_VF = abs(fft(ecg_VF.ecg(1:N))) .^ 2;
plot(((0:N-1)/N-0.5)*Fs_VF, fftshift(dsp_ecg_VF));
xlim([0;150]);
ylim([0;10^10]);

title('Power spectrum of ECG signal with VF');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot 223
N = Fs_SSS * s; % 4s
dsp_ecg_SSS = abs(fft(ecg_SSS.ecg(1:N))) .^ 2;
plot(((0:N-1)/N-0.5)*Fs_SSS, fftshift(dsp_ecg_SSS));
xlim([0;150]);
ylim([0;10^10]);

title('Power spectrum of ECG signal with SSS');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot 224
N = Fs_PVC * s; % 4s
dsp_ecg_PVC = abs(fft(ecg_PVC.ecg(1:N))) .^ 2;
plot(((0:N-1)/N-0.5)*Fs_PVC, fftshift(dsp_ecg_PVC));
xlim([0;150]);
ylim([0;10^10]);

title('Power spectrum of ECG signal with PVC');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

%% 4 Detection of P, QRS and T waves
%% 4.1 QRS
%% 4.1.1 R wave detection
clear all;

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

title('Time evolution of a normal ECG signal by the power of 3');
legend('ECG signal','R waves');
xlabel('Time (s)'); ylabel('Amplitude');

% ECG 2
subplot 312
N = Fs2 * s;
[pks,locs] = FindPeaks(ecg2P3(1:N), ratio, Fs2); % locate the local maxima

plot((0:N-1)/Fs2, ecg2P3(1:N));
hold on
plot((locs-1)/Fs2, pks,'+r');

title('Time evolution of a normal ECG signal by the power of 3');
legend('ECG signal','R waves');
xlabel('Time (s)'); ylabel('Amplitude');

% ECG 3
subplot 313
N = Fs3 * s;
[pks,locs] = FindPeaks(ecg3P3(1:N), ratio, Fs3); % locate the local maxima

plot((0:N-1)/Fs3, ecg3P3(1:N));
hold on
plot((locs-1)/Fs3, pks,'+r');

title('Time evolution of a normal ECG signal by the power of 3');
legend('ECG signal','R waves');
xlabel('Time (s)'); ylabel('Amplitude');

clear ecg1P3;
clear ecg2P3;
clear ecg3P3;
clear locs;
clear pks;
clear N;
clear S;
clear ratio;
%% Method of the derivation
s = 4; % seconds of samples
ratio = 70; % percent used to detect peaks

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
clear S;
clear locs;
clear pks;
clear ratio;
%% Pan and Tompkins algorithm
%% A)
clear all;
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
s = 10;
N = Fs1 * s; % S seconds of samples

figure(3);

subplot 231
plot((0:N-1)/N*s, ecg_n1.ecg(1:N));
title('Time evolution of ECG normal signal 1');
xlabel('Time (s)');
ylabel('Amplitude');

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
title('Time evolution of ECG normal signal 1 with low-pass filter');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 235
sp_y = abs(fft(y));
plot(((0:N-1)/N-0.5)*Fs1, fftshift(sp_y));
xlim([-50,50]);
title('Spectrum of ECG normal signal 1 with low-pass filter');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% High-pass filter
y1 = filter(B_hp, A_hp, y);

subplot 233
plot((0:N-1)/N*s, y1);
title('Time evolution of ECG normal signal 1 with band-pass filter');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 236
sp_y1 = abs(fft(y1));
plot(((0:N-1)/N-0.5)*Fs1, fftshift(sp_y1));
xlim([-50,50]);
title('Spectrum of ECG normal signal 1 with band-pass filter');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

clear sp_ecg_n1;
clear A_hp B_hp;
clear A_lp B_lp;

%% B)

B = [ 1 2 0 -2 -1 ] * Fs1 / 8;
A = [ 1 ];

figure(1);
freqz(B,A, 10000);
title('Filter amplitude and phase');

figure(2);

subplot 221
plot((0:N-1)/N*s, y1);
title('Time evolution of ECG normal signal 1 with band-pass filter');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 223
plot(((0:N-1)/N-0.5)*Fs1, fftshift(sp_y1));
xlim([-50,50]);
title('Spectrum of ECG normal signal 1 with band-pass filter');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

% Filter
y2 = filter(B, A, y1);

subplot 222
plot((0:N-1)/N*s, y2);
title('Time evolution of ECG normal signal 1 filtered');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 224
sp_y2 = abs(fft(y2));
plot(((0:N-1)/N-0.5)*Fs1, fftshift(sp_y2));
xlim([-50,50]);
title('Spectrum of ECG normal signal 1 filtered');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

clear A B;
%% C)

y3 = y2 .^ 2;

figure(1);

subplot 221
plot((0:N-1)/N*s, y2);
title('Time evolution of ECG normal signal 1 filtered');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 223
plot(((0:N-1)/N-0.5)*Fs1, fftshift(sp_y2));
xlim([-50,50]);
title('Spectrum of ECG normal signal 1 filtered');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot 222
plot((0:N-1)/N*s, y3);
title('Time evolution of ECG normal signal 1 filtered and sqared');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 224
sp_y3 = abs(fft(y3));
plot(((0:N-1)/N-0.5)*Fs1, fftshift(sp_y3));
xlim([-50,50]);
title('Spectrum of ECG normal signal 1 filtered');
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

clear N2 t;
%% E)

threshold = 4 * 10^11;

t = Smwi .* (Smwi >= threshold);
dt = diff(t);
D = [dt; 1:length(dt)];
figure(2);
plot(D(1,:));
select = D(1,:) >= threshold;
D(:,select==0) = [];

E = [-1*dt; 1:length(dt)];
select = E(1,:) >= threshold;
E(:,select==0) = [];
E(1,:) = -1* E(1,:);

figure(2);
plot(dt)
hold on
plot(D(2,:),D(1,:), '+r');
plot(E(2,:),E(1,:), '+g');
legend('Diff(ecg)','Maxima','Minima');
xlabel('Time (s)'); ylabel('Amplitude');

x = [D(2,:) E(2,:)]; % groups the locations
x = sort(x); % and sort them

if(dt(x(1)) < 0) % in case of the first data is not a maxima
    I=2:2:(length(x)-1);
else
    I=1:2:(length(x)-1);
end
R=[];
for i=I
    y = [ ecg_n1.ecg(x(i):x(i+1)); x(i):x(i+1) ]; % store in y the ecg data of the peak
    [maxi, id] = max(y, [], 2); % and get the peak and his location
    R = [ R [maxi(1); y(2,id(1))] ]; % store the results
end

peaks = R(1,:);
location = R(2,:);

figure(3);
plot(ecg_n1.ecg(1:N));
hold on;
plot(location, peaks, '+r');


%% 4.1.2 Q and S waves detection
clear all;
ecg_n1 = load('ecg/ecg_normal_1.mat');
Fs1 = ecg_n1.Fs;

s = 4;
ratio = 50;

N = Fs1 * s;
[pks,locs] = DerivMeth(ecg_n1.ecg(1:N), ratio, Fs1);

figure(1);
plot((0:N-1)/Fs1, ecg_n1.ecg(1:N));
hold on
plot((locs-1)/Fs1, pks,'+r');
title('Time evolution of a normal ECG signal');
legend('ECG signal','R waves');
xlabel('Time (s)'); ylabel('Amplitude');

%% QRS research

d = diff(ecg_n1.ecg(1:N));

M = 25; % Maximum of points to find the 0
Q = [];
R = locs;
S = [];

for i=1:length(R)
    % Q
    for j=1:M
        if(d(R(i) - j) <= 0)
            Q = [ Q (R(i)-j+1) ];
            break;
        end
    end
    % S
    for j=1:M
        if(d(R(i) + j) >= 0)
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
text(Q(1)/Fs1,ecg_n1.ecg(Q(1))*1.1,'Q'), text(R(1)/Fs1,ecg_n1.ecg(R(1))*1.03,'R'), text(S(1)/Fs1,ecg_n1.ecg(S(1))*1.1,'S')

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

[Q, R, S] = QRS(ecg_n1.ecg(1:N), ratio, ratio_deviation, Fs1);

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

[Q, R, S] = QRS(ecg_n2.ecg(1:N), ratio, ratio_deviation, Fs2);

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

[Q, R, S] = QRS(ecg_n3.ecg(1:N), ratio, ratio_deviation, Fs3);

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
%% 4.2.1 Filter
clear all;
ecg_n1 = load('ecg/ecg_normal_1.mat');
Fs1 = ecg_n1.Fs;

N = 1000;

% Differiator
B_g1 = [ 1 0 0 0 0 0 -1 ];
A_g1 = [ 1 ];

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

[pks,R] = DerivMeth(ecg_n1.ecg(1:N), ratio, Fs1);

R1 = R(2)
R2 = R(3)

E = R2 - R1
N = floor(E * 0.7)

data = ecg_n1.ecg(R1+Fs1/20:N+R1);

y = filter(B_g2, A_g2, data);
y = y(9:end);
y1 = filter(B_g1, A_g1, y);
y1 = y1(7:end);
figure(3);
subplot 211
plot(y1);
title('Time evolution of normal ECG signal 1 filtered');
xlabel('Time (s)');
ylabel('Amplitude');

subplot 212
sp = abs(fft(y1));
plot(fftshift(sp));
%xlim([-50,50]);
title('Spectrum of ECG normal signal 1 filtered');
xlabel('Frequency (Hz)');
ylabel('Amplitude')

%% Detection of the T waves
figure(1);

subplot 312
plot((0:length(y)-1),y);

subplot 313
plot((0:length(y1)-1),y1);
hold on

d = [y1, 1:length(y1)];
[M, idM] = max(d, [], 2);
[m, idm] = min(d, [], 2);
plot(idM-1, M, 'r+');
plot(idm-1, m, 'r+');

p = y1(idM:idm);
plot((idM-1:idm-1),p, 'r');

%find zero crossings
t1=p(1:end-1);
t2=p(2:end);
tt=t1.*t2;
indx=find(tt<0);

plot(indx+idM-2, p(indx),'+g');
N = length(data)

gd = 6;

subplot 311
plot((0:N-1),data);
hold on
plot(indx+idM +gd -1, data(indx+idM+gd), '+g');
%% P wave

figure(1);
subplot 411
    plot((0:999), ecg_n1.ecg(1:1000));
hold on

N = floor(E * 0.7)
data = ecg_n1.ecg(N+R1:R2-Fs1/20);

    plot((N-1+R1):(R2-1-Fs1/20), data, 'r');
    subplot 412
    plot((0:length(data)-1), data);
y = filter(B_g1, A_g1, data)
y = y(9:end);
    subplot 413
    plot((0:length(y)-1),y);
y1 = filter(B_g2, A_g2, y);
y1 = y1(7:end);
    subplot 414
    plot((0:(length(y1)-1)),y1);
d = [y1, 1:length(y1)];
[M, idM] = max(d, [], 2);
[m, idm] = min(d, [], 2);
p = y1(idM:idm);
    subplot 414

%find zero crossings
t1=p(1:end-1);
t2=p(2:end);
tt=t1.*t2;
indx=find(tt<0);
gd = 5;

figure(2);
plot((0:length(data)-1),data);
hold on
plot(indx+idM +gd -1, data(indx+idM+gd), '+g');


%% PQRST display of the 3 normal ECG signals
clear all;

ecg_n1 = load('ecg/ecg_normal_1.mat');
ecg_n2 = load('ecg/ecg_normal_2.mat');
ecg_n3 = load('ecg/ecg_normal_3.mat');

Fs1 = ecg_n1.Fs;
Fs2 = ecg_n2.Fs;
Fs3 = ecg_n3.Fs;

s = 50; % S seconds of samples

% ECG 1
subplot 311
N = Fs1 * s;

[ P, Q, R, S, T ] = PQRST(ecg_n1.ecg(1:N), Fs1);

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

[ P, Q, R, S, T ] = PQRST(ecg_n2.ecg(1:N), Fs2);

plot((0:N-1)/Fs2, ecg_n2.ecg(1:N));
hold on;
plot((P-1)/Fs2, ecg_n2.ecg(P),'+y');
plot((Q-1)/Fs2, ecg_n2.ecg(Q),'+g');
plot((R-1)/Fs2, ecg_n2.ecg(R),'+r');
plot((S-1)/Fs2, ecg_n2.ecg(S),'+k');
plot((T-1)/Fs2, ecg_n2.ecg(T),'+r');

title('Time evolution of normal ECG signal (1)');
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

[ P, Q, R, S, T ] = PQRST(ecg_n3.ecg(1+offset:N+offset), Fs3);

plot((0+offset:N-1+offset)/Fs3, ecg_n3.ecg(1+offset:N+offset));
hold on;
plot((P-1+offset)/Fs3, ecg_n3.ecg(P+offset),'+y');
plot((Q-1+offset)/Fs3, ecg_n3.ecg(Q+offset),'+g');
plot((R-1+offset)/Fs3, ecg_n3.ecg(R+offset),'+r');
plot((S-1+offset)/Fs3, ecg_n3.ecg(S+offset),'+k');
plot((T-1+offset)/Fs3, ecg_n3.ecg(T+offset),'+r');

title('Time evolution of normal ECG signal (1)');
xlabel('Time (s)');
ylabel('Amplitude');
text((P(1)+offset)/Fs3,ecg_n3.ecg(P(2)+offset)*-10,'P');
text((Q(1)+offset-40)/Fs3,ecg_n3.ecg(Q(2)+offset)*1.5,'Q');
text((R(1)+offset)/Fs3,ecg_n3.ecg(R(2)+offset)*1.1,'R');
text((S(1)+offset+10)/Fs3,ecg_n3.ecg(S(2)+offset)*0.8,'S');
text((T(1)+offset)/Fs3,ecg_n3.ecg(T(2)+offset)*1.5,'T');

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

N = 20 * Fs1;
CardiacRhythm(ecg_n3.ecg(1:N), Fs3)

% [ T1, B1 ] = FindPathologies(ecg_n1.ecg(1:N), Fs1)
% [ T2, B2 ] = FindPathologies(ecg_n2.ecg(1:N), Fs1)
% [ T3, B3 ] = FindPathologies(ecg_n3.ecg(1:N), Fs1)
% 
% [ T_AF, B_AF ] = FindPathologies(ecg_AF.ecg(1:N), Fs1)
% [ T_VF, B_VF ] = FindPathologies(ecg_VF.ecg(1:N), Fs1)
% [ T_SSS, B_SSS ] = FindPathologies(ecg_SSS.ecg(1:N), Fs1)
% [ T_PVC, B_PCV ] = FindPathologies(ecg_PVC.ecg(1:N), Fs1)


ratio = 40;
ratio_d = 1/12; % Fs * ratio_deviation : maximum deviation in point, to detect Q and S from R
% Q R S
N = 10 * Fs1;
[ Q, R, S ] = QRS(ecg_n1.ecg(1:N), ratio, ratio_d, Fs1);
N = length(R) - 1;
Ds = diff(R);
Dbpm = (60 * Fs1) * Ds.^-1;
subplot 211
plot((0:N-1), Ds, '+-');
subplot 212
plot((0:N-1), Dbpm, '+r-');

figure(2)s
% Q R S
N = 10 * Fs_SSS;
[ Q, R, S ] = QRS(ecg_SSS.ecg(1:N), 45, ratio_d, Fs_SSS);
N = length(R) - 1;w
Ds = diff(R);
Dbpm = (60 * Fs_SSS) * Ds.^-1;
subplot 211
plot((0:N-1), Ds, '+-');
subplot 212
plot((0:N-1), Dbpm, '+r-');












