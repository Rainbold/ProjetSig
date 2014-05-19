function [ n ] = CardiacRhythm(ecg, Fs)
   
    ratio = 40;
    ratio_d = 1/12; % Fs * ratio_deviation : maximum deviation in point, to detect Q and S from R
    % Q R S
    [ Q, R, S ] = QRS(2, ecg, Fs, ratio, ratio_d);
    
    N = length(R) - 1;
    
    D = (1/N) * sum(diff(R));
    
    n = 60 * Fs / D;
    
end