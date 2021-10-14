function Lpp = calculate_Lpp(time, Kp)
% Calcaulte plasma pause location according to Carpenter and Anderson [1992]

    Lpp = zeros(size(time));
    for it = 1:length(time)
        [tmp, first_idx] = min(abs(time - (time(it) - 1)));
        Kp24 = max(Kp([first_idx, it]));
        Lpp(it) = 5.6 - 0.46*Kp24;
    end
return