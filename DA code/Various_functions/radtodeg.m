%% degtorad
%
%
% Written by Adam Kellerman, Nov 2012. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Input
% 
% Name      Description                                                 
% -----------------------------------------------------------------------------------------
% angle         in radians

%% OUTPUT
% -----------------------------------------------------------------------------------------
% angle in degrees

%clear all
function [deg]=radtodeg(ang)
deg=ang*180./pi();
end
