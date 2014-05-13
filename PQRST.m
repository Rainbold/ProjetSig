function [ P, Q, R, S, T ] = PQRST( ecg, ratio, ratio_d, Fs )
    % Q R S
    [ Q, R, S ] = QRS(ecg, ratio, ratio_d, Fs);

    % Differiator
    B_g1 = [ 1 0 0 0 0 0 -1 ];
    A_g1 = [ 1 ];

    % Low-pass filter
    B_g2 = [ 1 0 0 0 0 0 0 0 -1 ];
    A_g2 = [ 1 -1 ];

    P = [];
    T = [];
    for i=1:length(R)-1
        R1 = R(i);
        R2 = R(i+1);

        E = R2 - R1;
        N = floor(E * 0.7);

        data = ecg(R1+Fs/20:N+R1);
        
        y = filter(B_g1, A_g1, data);
        y = y(7:end);
        y1 = filter(B_g2, A_g2, y)
    
        data2 = [y1, 1:length(y1)]
                figure(9)
        plot(data2);
        
        [M, idM] = max(data2, [], 2)
        [m, idm] = min(data2, [], 2)
        
        portion = y1(idM:idm)

        %find zero crossings
        t1=portion(1:end-1);
        t2=portion(2:end);
        tt=t1.*t2;
        indx=find(tt<0);
        
        T = [ T (indx+idM+R1)];
    end
    
end

