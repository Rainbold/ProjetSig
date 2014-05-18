function [ Tachy, Brady, Ectopic ] = FindPathologies( ecg, Fs, ratio)

    n = CardiacRhythm(ecg, Fs);
    [~,R] = DerivMeth(ecg_n2.ecg(1:N), ratio, Fs2);
    
    % Tachycardia
    if(n > 100)
       Tachy = true; 
    else
        Tachy = false;
    end
    
    % Bradicardia
    if(n < 60)
        Brady = true;
    else
        Brady = false;
    end
    
    % Ectopic
    dr = diff(R);
    ddr = diff(dr);
    threshold = 15;
    if(ddr(ddr > threshold) ~= 0)
       Ectopic = true; 
    else
        Ectopic = false;
    end
    
end

