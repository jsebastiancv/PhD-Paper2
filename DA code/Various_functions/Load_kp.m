%% Load_kp

% Outputs: kp_time, kp_index

% loads Kp from Jan 1932 to the most recent update.

% kp_time - matlab serial date number
% kp_index - Kp index in x10 format (OMNI)
% Satellite_data_dir - path to the mounted and marked /data/ directory on the data server (yssvnl)

% The Kp times here are centered on the 3-hr interval of interest !!!!!!
% this is done so that interpolation routines using nearest neighbour will
% find the correct Kp index for the period of interest.

%% this Kp file is made from a script:
% yssvnl.igpp.ucla.edu:/data/Kp/Kp_archive/Process_Kp

% Written by Adam Kellerman, Jan 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Work is supported by:
% LANL grant 12-LR-235337
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kpfile = [Satellite_data_dir,'Kp/Kp_archive/kp_all.mat'];
load(kpfile)
fprintf('Note that Kp times are centered on each 3-hr interval\n');

return

