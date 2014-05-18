function [ StdDev ] = HRV( ecg, Fs, ratio, ratio_d )

    [Q, R, S] = QRS(ecg, Fs, ratio, ratio_d);
    
    dr = diff(R);
    
    % Standard deviation
    StdDev = sqrt( (sum(( dr - mean(dr) ).^2)) / length(dr) );
    
    % Histogram of occurences
    figure(1)
    hist(dr);
    title('Histogram of occurences');
    
    % Display of process v(t)
    N = length(dr)
    v = dr(1:N-1) + ( (dr(2:N) - dr(1:N-1)) .* ((2:N) - (1:N-1)).^-1);
    figure(2);
    plot(v);
    title('Process v(t)');
    

end

