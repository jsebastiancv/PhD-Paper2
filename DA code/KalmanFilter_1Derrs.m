%% KalmanFilter - analysis phase
% No parameterization
%
%
% Written by Marianne Daae Oct 14, 2009
% slightly modified for pll compatibility by A. Kellerman October, 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Explaining Variables
%
% Name      Description                                                 Dim
% -----------------------------------------------------------------------------------------
% nx        spatial dimension                                           1
% it        current time index                                          1
% nobs      Number of observations available for this analysis          1
% Me        Model error                                                 1
% Re        Observation error                                           1
% Xf        Forecast                                                    nx
% d         Innovation vector                                           nx
% Kd        Weighted innovation                                         nx
% Xo        Observations                                                nobs
% Q         Model error covariance matrix                               nx x nx
% Pf        Forecast error covariance matrix                            nx x nx
% M         Model Matrix Operator                                       nx x nx
% K         Kalman matrix, or optimal gain matrix                       nx x nobs
% H'        Model operator, transforms observations into model space    nx x nobs
% H         Model operator, transforms forecast into observation space  nobs x nx
% R         Observation error covariance matrix                         nobs x nobs


% OUTPUT
% -----------------------------------------------------------------------------------------
% Xa        Analysis                                                    nx
% Pa        Analysis error covariance matrix                            nx x nx


%function [Xa,Pf,Q,R,d,Kd]=KalmanFilter_1Dnew(Xf,Me,M,Pa,SatDat,rr)
function [Xa,Pf]=KalmanFilter_1Derrs(Xf,Me,M,Pa,SatDat,rr)

% global no_data_asim

[Xo,Re,H,nobs] = get_obs2D(SatDat,rr);

fmin = 1.e-5;
if(size(Xf,1)==1)
    Xf=Xf';
end
nx=length(Xf);
q=(Me.*Xf).^2;
q(q <= fmin)=fmin;
R = NaN;
Q  = diag(q);  %model error is fraction of value so this is ok
%Q = q*q';
%km1=diag(Q,-1);
%k0=diag(Q);
%kp1=diag(Q,1);
%Q = diag(km1,-1)+diag(k0)+diag(kp1,1);
Pf = M*Pa*M'+Q;

if nobs>0
    r = (Re.*Xo).^2;
    R  = diag(r); %make covariance
    rc_R = rcond(R); %check condition
     if rc_R > 0.001
        
        %R  = diag(Re);
        K  = Pf*H'/(H*Pf*H'+R);
        %K'
        %fprintf('Rcond: %s \n',num2str(rc_R));
        %fprintf('nobs: %s \n',num2str(nobs));
        d  = Xo - H*Xf;
        Hd = H'*d;
        Od = H'*Xo;
        Kd = K*d;
        Xa = Kd+Xf;
        Pf = (eye(nx)-K*H)*Pf;
    
     else
       %  fprintf('Bad condition on R: %s\n',num2str(rc_R));
         Xa = Xf;
         d = [];
         Kd = [];
     end

else
    Xa = Xf;
    d = [];
    Kd = [];
end

return
