function [ StdDev, dr, v, bpm ] = HRV(method, ecg, Fs, ratio, ratio_d )

    [~, R, ~] = QRS(method, ecg, Fs, ratio, ratio_d);
    
    dr = diff(R);
    
    % Standard deviation
    StdDev = sqrt( (sum(( dr - mean(dr) ).^2)) / length(dr) );
    
    % Display of process v(t)
    N = length(dr);
    v = [];
    dr = diff(R);
    for i=1:N-1
        
        for j=R(i):R(i+1)
            v = [ v (dr(i) + ( (dr(i+1)-dr(i) ) .* (R(i+1) - R(i)).^-1 ) * (j-R(i))) ];
        end
    end
    N = 10^5;
    sp = fftshift(abs(fft((v-mean(v)),N)));
    half = sp((N/2+1):end);
    % N/2 -> Fs/2 Hz
    % X -> 0.4Hz
    % X = N/2 * 0.4 / (Fs/2)

    data = half(1:((N/2) * 0.4 / (Fs/2)));
    d = diff(data);
    dd = diff(d);
    t1=d(1:end-1);
    t2=d(2:end);
    tt=t1.*t2;
    indx=find(tt<0);
    indx(dd(indx) > 0) = [];
    
    if(length(indx) >= 4)
        bpm = indx(4) * 60 * Fs / N;
    else
        bpm = 0;
    end
end