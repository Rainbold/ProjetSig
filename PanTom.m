function [ R ] = PanTom( ecg, Fs, ratio )

    % Low-pass filter
    B_lp = zeros(1,13);
    B_lp(1) = 1;
    B_lp(7) = -2;
    B_lp(13) = 1;

    A_lp = [ 1 -2 1 ];

    % High-pass filter
    B_hp = zeros(1,33);
    B_hp(1) = -1;
    B_hp(17) = 32;
    B_hp(18) = -32;
    B_hp(33) = 1;

    A_hp = [ 1 1 ];

    % Low-pass filter
    y = filter(B_lp, A_lp, ecg);
    
    % High-pass filter
    y1 = filter(B_hp, A_hp, y);

    B = [ 1 2 0 -2 -1 ] * Fs / 8;
    A = 1;

    % Filter
    y2 = filter(B, A, y1);
    
    y3 = y2 .^ 2;
    
    N2 = 2/10 * Fs; % Value for the length of the window
    Smwi = [];

    for n=N2:length(y3);
        t = (1/N2) * sum(y3((n-N2+1):n));
        Smwi = [ Smwi t ];
    end
    
    max_s = max(Smwi); % Maximum of integrated ECG

    [pks, maxima] = FindPeaks(Smwi, ratio, Fs);

    % Detection of the beginin of each peaks
    low_threshold = 0.05 * max_s;
    data = Smwi;
    data(data > low_threshold) = 0;
    ds = diff(data);

    [pks2, minima] = FindPeaks(ds, 10, Fs);

    minima = minima(1:2:end);
    pks2 = pks2(1:2:end);

    % Detection of R peaks
    x = [minima maxima]; % groups the minmas/maximas
    x = sort(x); % and sort them

    R=[];
    if(numel(x) ~= 0) % Protection useful for the GUI
        if(x(1) ~= minima(1)) % in case of the first data is not a minima
            I=2:2:(length(x)-1);
        else
           I=1:2:(length(x)-1); % used to get every couple maxima-minima
        end

        delay = N2/2;

        for i=I
            y = [ ecg(x(i)+delay:x(i+1)+delay); x(i)+delay:x(i+1)+delay ]; % store in y the ecg data of the peak and their #
            [maxi, id] = max(y, [], 2); % get the peak and his location with their #
            R = [ R [maxi(1); y(2,id(1))] ]; % store the results
        end
    end
    if(size(R) ~= [0,0])
        R = R(2,:);
    else
        R = [];
    end
end