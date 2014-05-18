function [ Q, R, S ] = QRS(method, ecg, Fs, ratio, ratio_d)
    Q = [];
    S = [];
    max_d = Fs * ratio_d;
    
    % R
    switch method
        case 1
            [pks,R] = FindPeaks(ecg.^3, ratio, Fs);
        case 2
            [pks,R] = DerivMeth(ecg, ratio, Fs);
        case 3
            [ R ] = PanTom(ecg, Fs, ratio);
    end
    
    d = diff(ecg);
    for i=1:length(R)
        % Q
        for j=1:max_d
            if((R(i) - j) > 0)
                if(d(R(i) - j) < 0)
                    Q = [ Q (R(i)-j+1) ];
                    break;
                end
            end
        end
        % S
        for j=1:max_d
            if((R(i)+j) <= length(d))
                if(d(R(i) + j) > 0)
                    S = [ S (R(i)+j) ];
                    break;
                end
            end
        end
    end
end

