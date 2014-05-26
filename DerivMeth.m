function [peaks, location] = DerivMeth(ecg, ratio, Fs)
    
    % Because noiseBL ecg data is in line and not in collun like everyelse
    % signals, we have to do this test.
    s = size(ecg);
    if(s(1) ~= 1)
       ecg = (ecg)'; 
    end

    ecg_d = diff(ecg); % derivative of ecg

    [pksM,locsM] = FindPeaks(ecg_d, ratio, Fs); % locate the local maximas
    [pksm,locsm] = FindPeaks(-1*ecg_d, ratio, Fs); % locate the local minimas
    pksm = -1 * pksm; % put back the minimas under zero

%debuging :
%     figure(9);
%     subplot 211
%     plot((1:length(ecg_d))/Fs, ecg_d);
%     hold on;
%     plot(locsM/Fs, pksM,'+r');
%     plot(locsm/Fs, pksm,'+g');
%     title('Time evolution of the differential of an ECG signal');
%     legend('Diff(ecg)','Maxima','Minima');
%     xlabel('Time (s)'); ylabel('Amplitude');
%     figure(1);

    x = [locsm locsM]; % groups the locations
    x = sort(x); % and sort them
    R=[];
    if(numel(x) ~= 0) % Protection useful for the GUI
        if(ecg_d(x(1)) < 0) % in case of the first data is not a maxima
            I=2:2:(length(x)-1);
        else
           I=1:2:(length(x)-1); % used to get every couple maxima-minima
        end
        
        for i=I
            y = [ ecg(x(i):x(i+1)); x(i):x(i+1) ]; % store in y the ecg data of the peak and their #
            [maxi, id] = max(y, [], 2); % get the peak and his location with their #
            R = [ R [maxi(1); y(2,id(1))] ]; % store the results
        end
    end
    if(size(R) ~= [0,0]) % Protection useful for the GUI
        peaks = R(1,:);
        location = R(2,:);
    else
        peaks = [];
        location = [];
    end
end

