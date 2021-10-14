function y=Kfunc(pc)
mc2 = 0.511;
y = ( sqrt( 1.0 + ( pc ./ mc2).^2 ) - 1) .* mc2;
return;
