function res = Y(alpha)
    T0 = 1.3802;
    T1 = 0.7405;
    res = 2.*(1.0 - sin(alpha)).*T0 + (T0 - T1).*(sin(alpha).*log(sin(alpha)) + 2.*sin(alpha) - 2.*sqrt(sin(alpha)));
return 
