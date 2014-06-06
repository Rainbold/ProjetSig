function [ AF, VF ] = AFVF( ecg, R, Fs )

    % AF
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
    if(max(abs(gamma)) >= 1000)
        AF = true;
    else
        AF = false;
    end
    
    % VF
    N = length(ecg);
    sp_ecg = abs(fft(ecg));
    sp_ecg(1) = 0;
    sp_ecg = fftshift(sp_ecg);

    N = length(sp_ecg);
    sp = sp_ecg(floor(N/2):end);
    id_max_sp = find(sp == max(sp));

    rhythm = id_max_sp * Fs / N * 60;
    
    if(rhythm >= 300)
        VF = true;
    else
        VF = false;
    end
end

