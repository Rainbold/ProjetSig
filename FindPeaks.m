function [peaks, location] = FindPeaks(ecg, ratio, Fs)

    max_peak = max(ecg);
    threshold = max_peak * ratio/100;

    D = [ecg; 1:length(ecg)]; % ecg data with # of samples
    select = D(1,:) >= threshold; % all data under threshold go to zero
    D(:,select==0) = []; % we remove the zeros

    P = [diff(D(2,:)); D(2,1:end-1)]; % differential of the # of samples to detect the different peaks
    R = [];
    limit = Fs / 10; % Limit to detect the different peaks in the derivation of the #
    d = 1;
    for j1=1:(length(D(2,:))-1)
        if(P(1,j1) >= limit || j1 == (length(D(2,:))-1)) % if the diff is more than 50 or if we are at the end of the diff data
            f = j1;
            L = [ D(1,d:f); D(2,d:f) ]; % ecg data of a peak
            [maxi, i] = max(L, [], 2); % get the value of the peak and his #
            R = [ R [maxi(1); L(2,i(1))] ]; % save results
            d = j1 + 1; % set d for the next peak
        end
    end
    if(size(R) ~= [0,0])
        peaks = R(1,:);
        location = R(2,:);
    else
        peaks = [];
        location = [];
    end
end