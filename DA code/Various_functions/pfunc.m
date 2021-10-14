function y=pfunc(K)
mc2 = 0.511;
y = sqrt( (K ./ mc2 + 1).^2 - 1) .* mc2;
return;
