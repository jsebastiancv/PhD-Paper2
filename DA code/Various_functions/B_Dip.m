function B_mag = B_Dip(R, lambda, B0)
    %%
    % B_DIP calculates dipole magnetic field in equatorial plane
    %
    % Input:
    %   R - radial distance, Re
    %   lambda - latitude, rad (lambda = 0 if not specified)
    %   B0 - (optional) magnetic field at Earth's surface in Gauss
    %
    % Output:
    %   B_mag - magnetic field in Gauss

    % Author: Dmitriy Subbotin
    % Email: subbotin@ucla.edu
    % Last change: 20014-01-01 by Dmitriy Subbotin
    % The work was supported by:
    %   LANL grant 12-LR-235337, PI Yuri Shprits
    %   NASA grant NNX09AF51G, PI Yuri Shprits
    
    if nargin < 2
        lambda = 0;
    end
    if nargin < 3
        B0 = 0.31;%2; %Gauss
    end

    %B_r = -2 * B0 ./ (L.^3) .* sin(lambda);
    %B_t = -B0 ./ (L.^3) .* cos(lambda);
    B_mag = B0 ./ (R.^3) .* sqrt(1 + 3 * sin(lambda).^2) ./ cos(lambda).^6;
    
end

