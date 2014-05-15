function [ Tachy, Brady ] = FindPathologies( ecg, Fs )

    n = CardiacRhythm(ecg, Fs);
    if(n > 100)
       Tachy = true; 
    else
        Tachy = false;
    end
    if(n < 60)
        Brady = true;
    else
        Brady = false;
    end
end

