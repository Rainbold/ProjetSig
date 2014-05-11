function [ Q, R, S ] = QRS(ecg, ratio, max_d)
    Q = [];
    R = [];
    S = [];

    % R
    [pks,R] = DerivMeth(ecg, ratio);

    d = diff(ecg);
    for i=1:length(R)
        % Q
        for j=1:max_d
            if(d(R(i) - j) < 0)
                Q = [ Q (R(i)-j+1) ];
                break;
            end
        end
        % S
        for j=1:max_d
            if(d(R(i) + j) > 0)
                S = [ S (R(i)+j) ];
                break;
            end
        end
    end
end

