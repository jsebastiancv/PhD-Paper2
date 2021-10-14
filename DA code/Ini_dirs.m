% Set up directory and variable structure

% Written by Adam Kellerman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Work is supported by:
% LANL grant 12-LR-235337

global DA_dir Satellite_data_dir VERB_dir IRBEM_dir

%DA_dir=[getenv('HOME'),'/projects/3d_verb_data_assimilation/']; % mark it!
DA_dir = [getenv('HOME'),'/Documents/VERB/22. 3D VERB DA_Spectrum_LT/']; 

% path to sshfs mounted data server root (yssvnl.epss.ucla.edu:/data/)
%Satellite_data_dir=([getenv('HOME'),'/data/']); % mark it! 
Satellite_data_dir = ['/export/.automounter/mag5/data/rbm-data/'];
%Satellite_data_dir = ('Y:'); 

% path to VERB code 2.0 root directory
VERB_dir = [getenv('HOME'),'/Documents/VERB/VERB 2.x/VERB open mp/']; % mark it!
%VERB_dir = ('C:\Users\jscv\Documents\MATLAB\VERB 2.0'); 

% path to MATLAB routines, with compiled IRBEM-4.4.0 code
%IRBEM_dir = (['/var/IRBEM-4.4.0/matlab/']); % mark it!
%IRBEM_dir = ('C:\Users\jscv\Documents\MATLAB\ONERA\matlab');
IRBEM_dir = [getenv('HOME'),'/Documents/VERB/ONERA/matlab/'];

    % you can of course change this given a new version, though some code changes may be required depending
    % on any updates

fprintf('Home directory set to %s\n',DA_dir);
fprintf('Data direcotry set to %s\n',Satellite_data_dir);
fprintf('VERB dir set to %s\n', VERB_dir);
fprintf('IRBEM dir set to %s\n', IRBEM_dir);

% check for existance. 

check_1 = exist(DA_dir,'dir');
check_2 = exist([Satellite_data_dir,'GOES/'],'dir'); % if sshfs has mounted then this should exist
check_3 = exist(VERB_dir,'dir');
check_4 = exist(IRBEM_dir,'dir');

if ~check_1 == 7 %if DA_dir is correct
    fprintf('Set home path: %s\n. Detected current working directory: %s\n',DA_dir,pwd);
    error('Your home directory is not correct, you must specify this in Ini_dirs.m');
end
if ~check_2 == 7 %if Data exists
    error('It does not look like you mounted yssvnl.epss.ucla.edu:/data correctly, or you path to your local repository is incorrect');
end
if ~check_3 == 7 %if VERB exists
    error('The path to the VERB code does not exist');
end
if ~check_4 == 7 %if IRBEM exists
    error('The path to IRBEM does not exist');
end

% add paths to MATLAB workspace
addpath([DA_dir,'Code/']);
addpath([DA_dir,'Code/Various_functions/']);
addpath(IRBEM_dir);
        