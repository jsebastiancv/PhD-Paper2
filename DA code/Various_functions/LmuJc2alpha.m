function alpha1 = LmuJc2alpha(L, mu, Jc)
%LMUJC2ALPHA convert L, mu and J to alpha
    [L_mesh, mu_mesh, Jc_mesh] = meshgrid(L, mu, Jc);

    mc2 = 0.511;
    for i = 1:length(L_mesh(:)) 
        alpha1(i) = fzero(@(alpha) (Y(alpha)./sin(alpha) - Jc_mesh(i) .* sqrt(L_mesh(i)) ./ sqrt(8.0 .* B0(1) .* mc2 .* mu_mesh(i))), [0.001, pi/2-0.001]);
    end
    alpha1 = squeeze(reshape(alpha1, size(L_mesh)));
    
return
