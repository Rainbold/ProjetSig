function [ n ] = CardiacRhythm(ecg, Fs)
   
    ratio = 40;
    ratio_d = 1/12; % Fs * ratio_deviation : maximum deviation in point, to detect Q and S from R
    % Q R S
    [ Q, R, S ] = QRS(ecg, ratio, ratio_d, Fs);
    
    N = length(R) - 1;
    
    D = (1/N) * sum(diff(R));
    
    n = 60 * Fs / D;
    
end