function [ P, Q, R, S, T ] = PQRST( method, ecg, Fs, ratio, ratio_d )
   
    % Because noiseBL ecg data is in line and not in column like everyelse
    % signals, we have to do this test.
    s = size(ecg);
    if(s(1) ~= 1)
       ecg = (ecg)'; 
    end

    % Q R S
    [ Q, R, S ] = QRS(method, ecg, Fs, ratio, ratio_d);

    % Differiator
    B_g1 = [ 1 0 0 0 0 0 0 0 -1 ];
    A_g1 = [ 1 -1 ];

    % Low-pass filter
    B_g2 = [ 1 0 0 0 0 0 -1 ];
    A_g2 = [ 1 ];

    P = [];
    T = [];
    for i=1:length(R)-1    
        R1 = R(i);
        R2 = R(i+1);

        E = R2 - R1;

        % P
        N = floor(E * 0.7);
        data = ecg(R1+N:R2-floor(Fs/20));
        y = filter(B_g1, A_g1, data);
        y = y(9:end);
        y1 = filter(B_g2, A_g2, y);
        y1 = y1(7:end);
        d = [y1; 1:length(y1)];
        [M, idM] = max(d, [], 2);
        [m, idm] = min(d, [], 2);
        if(numel(idM)~=0 && numel(idm)~=0)
            p = y1(idM(1):idm(1));
            %find zero crossings
            t1=p(1:end-1);
            t2=p(2:end);
            tt=t1.*t2;
            indx=find(tt<0);
            if(numel(indx) > 0)
                indx = indx(1);
            end
            group_delay = 6;
            P = [ P (indx + R1+N + idM(1) + group_delay)];
        end
        
        % T
        N = floor(E * 0.7);
        data = ecg(R1+floor(Fs/20):N+R1);
        y = filter(B_g1, A_g1, data);
        y = y(9:end);
        y1 = filter(B_g2, A_g2, y);
        y1 = y1(7:end);
        
        d = [y1; 1:length(y1)];
        [M, idM] = max(d, [], 2);
        [m, idm] = min(d, [], 2);
        if(numel(idM)~=0 && numel(idm)~=0)
            p = y1(idM(1):idm(1));
            %find zero crossings
            t1=p(1:end-1);
            t2=p(2:end);
            tt=t1.*t2;
            indx=find(tt<0);
            if(numel(indx) > 0)
                indx = indx(1);
            end
            group_delay = 6;
            T = [ T (indx + R1+floor(Fs/20) + idM(1) + group_delay)];
        end
    end
    
end

