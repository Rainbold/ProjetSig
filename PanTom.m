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

    threshold = 4 * 10^11 * 2 * (ratio/100);

    t = Smwi .* (Smwi >= threshold);
    dt = diff(t);
    D = [dt; 1:length(dt)];
    select = D(1,:) >= threshold;
    D(:,select==0) = [];

    E = [-1*dt; 1:length(dt)];
    select = E(1,:) >= threshold;
    E(:,select==0) = [];
    E(1,:) = -1* E(1,:);

    x = [D(2,:) E(2,:)]; % groups the locations
    x = sort(x); % and sort them

    R=[];
    if(numel(x) ~= 0)
        if(dt(x(1)) < 0) % in case of the first data is not a maxima
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
        R = R(2,:);
    else
        R = [];
    end
end