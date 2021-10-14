function [meanval]=wmean(input,weights)

zz= find(~isnan(input));
input = input(zz);
weights = weights(zz);

nel = length(input);

num = sum(weights.*input);
den = sum(weights);

meanval = num./den;
