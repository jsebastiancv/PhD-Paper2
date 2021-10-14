function PSD_SS = steady_state(L1d, D1d, tau)
    %%
    % STEADY_STATE solves d/dL DLL/L^2 df/dL + f/tau= 0 equation
    %
    % PSD_SS = steady_state(L1d, D1d, tau)
    %
    % Input:
    %   L1d - 1d L-star array
    %   D1d - 1d radial diffusion coefficient
    %   tau - electron lifetime
    %
    % Output:
    %   PSD_SS - computed 1d PSD

    % Author: Yuri Shprits
    % Email: shprits@mit.edu
    % Last change: 20014-01-01 by Dmitriy Subbotin
    % The work was supported by:
    %   LANL grant 12-LR-235337, PI Yuri Shprits
    %   NASA grant NNX09AF51G, PI Yuri Shprits
    
    Am1 = zeros(size(L1d));
    A0 = zeros(size(L1d));
    Ap1 = zeros(size(L1d));
    
    if (length(tau) == 1)
        tau = repmat(tau, size(L1d));
    end
    
    for iL = 2:length(L1d)-1
        Gdh = 1.0/((L1d(iL+1) - L1d(iL-1))/2) .* (L1d(iL).^2);
        dfdL_L = ((D1d(iL) + D1d(iL-1))/2)/(L1d(iL) - L1d(iL-1))/((((L1d(iL) + L1d(iL-1))/2).^2));
        dfdL_R = ((D1d(iL+1) + D1d(iL))/2)/(L1d(iL+1) - L1d(iL))/((((L1d(iL+1) + L1d(iL))/2).^2));
        Am1(iL) = Gdh * dfdL_L;
        A0(iL) = -Gdh * (dfdL_L + dfdL_R) - 1./tau(iL);
        Ap1(iL) = Gdh * dfdL_R;
    end
    RHS = zeros(size(L1d));
    RHS(1) = 1e-99;
    RHS(end) = 1;
    
    A0(1) = 1;
    A0(end) = 1;

    A = diag(A0) + diag(Am1(2:end), -1) + diag(Ap1(1:end-1), 1);
    
    PSD_SS = A\RHS(:);
    %PSD_SS = inv(A) * RHS(:);

end
