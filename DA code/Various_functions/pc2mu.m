function y = pc2mu (L, pc, alpha)
   mc2 = 0.511;
   y = pc.^2 .* sin(alpha).^2 ./ (B0(L) * 2 * mc2);
return;
