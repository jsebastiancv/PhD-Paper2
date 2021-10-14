% This is the main script from which to run the Reanalysis Code.
% the purpose of this script is to create reanalysis files for monthly periods
%
% Set up the Ini_dirs.m file first, which covers the location of the processed
% spacecraft dataset and several other important variables related to folder structure
%
% Specify the reanalysis parameters in the setup_DA3D.
% Please read through each option carefully before running. 

% Written by Adam Kellerman. 2012-2016.
% Last change 2016-01-21

% Work is supported by:
% LANL grant 12-LR-235337

clear all

setenv('LD_LIBRARY_PATH', getenv('PATH'))

%% initialize directories
Ini_dirs

%% build array with upper L boundary condition derived from GOES data

build_spectrum_from_data_with_logs_missing_month_5th_test

%% run setup code
setup_DA3D

% create simulation directory if necessary
make_sim_dir
%% choose time period, files are created monthly, starting with a SS initial condition for the very first time point, then last psd thereafter

% These are the start and end months, inclusive, for which to run the reanalysis
styr = 2015;
stmo = 09;
% stmo = 10;
enyr = 2016;
enmo = 09;

if save_verb_forecast
    if load_errors_from_file
        ini_model_err = false;
    else
        ini_model_err = true;
    end
    if timestep ~= 1
        error('The default time step for error determination is one hour, this is hard coded atm')
    end
else
    ini_model_err = false;
end
%% if not starting from steady state then load the restart_file
if ~use_SS
    fprintf('loading restart file: %s...',restart_file);
    load(restart_file);
    if use_restart_time
        [~,inx] = min(abs(SimTime - restart_time));
        last_PSD = squeeze(SimPSD(inx,:,:,:));
    else
        last_PSD = squeeze(SimPSD(end,:,:,:));
    end
    fprintf('done\n');
end

for iyear = styr:enyr
    if iyear == styr
        smo = stmo;
    else
        smo = 1;
    end
    if iyear == enyr
        emo = enmo;
    else
        emo = 12;
    end
    for imo = smo:emo
        dim = find_DIM(iyear, imo); % days in month
        start_time = datenum(iyear,imo,1); % inclusive, initial condition
        end_time = datenum(iyear,imo,dim)+1; % inclusive, as we want an initial condition beginning
                                                % on the start date for the next loop
        yrmo = datestr(start_time,'yyyymm');

        DA3D_Neumann_data_from_previous_run
    end
end
