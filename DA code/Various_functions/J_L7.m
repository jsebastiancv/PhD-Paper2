function f=J_L7(K, x1, y1, x2, y2)
%   x1=0.2;
%   y1=2.e3;
%  x2=1.;
%   y2=7.;
   bcoef = log(y1/y2)/(x2-x1);
   acoef = y2*exp(bcoef*x2);
   f=acoef*exp(-bcoef*K);
return

