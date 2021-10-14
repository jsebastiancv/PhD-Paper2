function res=J_L7_corrected(K)
x = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 2e-1, 1e0, 3e0, 1e1];%, 2e1];
y = [3e9,  2e7,  9e6,  3e5,  1e4,  2e3, 7e0, 3e-3, 1e-3];%, 1e-3];

i=2;
res = [];
for each_k = K
    while ((x(i) < each_k) && (i < length(x))) i = i + 1; end;
    f = J_L7(each_k, x(i-1), y(i-1), x(i), y(i));
    res=[res, f];
end

return 