function [peaks, location] = DerivMeth(ecg, ratio, Fs)
    
    ecg_d = diff(ecg);

    [pksM,locsM] = FindPeaks(ecg_d, ratio, Fs); % locate the local maximas
    [pksm,locsm] = FindPeaks(-1*ecg_d, ratio, Fs); % locate the local minimas
    pksm = -1 * pksm;

%debuging :
%     figure(9);
%     
%     plot(ecg_d);
%     hold on;
%     plot(locsM, pksM,'+r');
%     plot(locsm, pksm,'+g');
%     title('Time evolution of the differential of an ECG signal');
%     legend('Diff(ecg)','Maxima','Minima');
%     xlabel('Time (s)'); ylabel('Amplitude');
%     figure(1);

    x = [locsm locsM]; % groups the locations
    x = sort(x); % and sort them
    R=[];
    if(numel(x) ~= 0)
        if(ecg_d(x(1)) < 0) % in case of the first data is not a maxima
            I=2:2:(length(x)-1);
        else
           I=1:2:(length(x)-1);
        end
        
        for i=I
            y = [ ecg(x(i):x(i+1)); x(i):x(i+1) ]; % store in y the ecg data of the peak
            [maxi, id] = max(y, [], 2); % and get the peak and his location
            R = [ R [maxi(1); y(2,id(1))] ]; % store the results
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

