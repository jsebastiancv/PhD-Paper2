function [DLL] = RadialDiffCoeff(Kp,L)
%RadialDiffCoeff - This function returns the radial diffusion coefficients
%for L between 2 and 8 and for a provided value of Kp. 
%   This function enables one to reproduce the first row of figures in this paper

DLL = 10^(0.506*Kp - 9.325)*(L.^10);

%To see plot, uncomment following lines.

%L = linspace(2,8);
%plot(L,1./DLL(2,L),'.-k'); 
%set(gca,'yscale','log','ylim',[1e-2 1e6],'ylabel','1/Dll','xlabel,'L');
end

