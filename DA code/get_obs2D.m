%% Finding irregular grided obsevations 
% 
%
% Written by Marianne Daae Oct 14, 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Explaining Variables
% 
% Name      Description                                                 Dim
% -----------------------------------------------------------------------------------------
% nx        spatial dimension                                           1
% nit       temporal dimenstion                                         1
% it        current time index                                          1
% nobs      Number of observations available for this analysis          1
% rr        Errors of Observation                                       5
% indobs    position of ovservations                                    nobs
% SatDat    All satelite data                                           structure
% tmpObs    Current satellite date                                      nit x nx


% OUTPUT
% -----------------------------------------------------------------------------------------
% Re        Observation error                                           1
% Xo        Observations                                                nobs
% H         Model operator, transforms forecast into observation space  nobs x nx 

%%
function [Xo,Re,H,nobs] = get_obs2D(SatDat,rr)


% global rr


tmpObs= reshape(SatDat',[],1);
nx=length(tmpObs);
indobs = find(~isnan(tmpObs));  
nobs   = length(indobs);                           
H      = zeros(nobs,nx);                          
Xo          = tmpObs(indobs);
Re = rr(indobs);
%if ~isempty(indobs)
%	Xo'
%	Re'
%end
for ik=1:nobs
  H(ik,indobs(ik))=1.; % hard wired at this point in time. We may in future start with fluxes and use this 
                        % to go to PSD
end

return
