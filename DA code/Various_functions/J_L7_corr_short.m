function res=J_L7_corr_short(K)

x_arr = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1];

res = [];
for X = K'
    if (X <= x_arr(2)) 
        f = 5.23491e9 * exp(-5.56737e4 * X); 
    elseif (X > x_arr(2) && X <= x_arr(3)) 
        f = 2.18556e7 * exp(-8.87231e2 * X); 
    elseif (X > x_arr(3) && X <= x_arr(4)) 
        f = 1.31331e7 * exp(-3.77911e2 * X); 
    elseif (X > x_arr(4) && X <= x_arr(5)) 
        f = 4.37770e5 * exp(-3.77911e1 * X); 
    elseif (X > x_arr(5) && X <= x_arr(6)) 
        f = 2.24153e4 * exp(-8.07159e0 * X); 
    elseif (X > x_arr(6)) 
        f = 1.87211e1 * exp(-9.83741e-1 * X); 
    else
        'ERROR';
        return; 
    end;
    res=[res, f];
end

return 