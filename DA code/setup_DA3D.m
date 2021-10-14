% Reanalysis Code. This code operates in a month-to-month format.
% Set up the Ini_dirs.m file first, which covers the location of the processed spacecraft dataset and 
% several other important variables
%
% Specify the start month, the grid format, the spacecraft to assimilate, and other options.
% Please read through each option carefully before running. 

% Written by Adam Kellerman. 2012-2016.
% Last change 2016-01-21

% Work is supported by:
% LANL grant 12-LR-235337

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin code options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time step for model and assimilation given in hours
    timestep = 1; %hours, this controls the assimilation and matrices output from VERB

%update step for reanalysis
    ndastep = 1; % iteration index
    % Note that by default this includes data in between the update times. If you want to just
    %assimilate data that matches the last time step of the simulation then specify here
        use_data_snapshot = false; %[false], true
        % true, will only assimilate data matching the last time step of the simulation t to t+1;

% number of cpu's to use if you have the parallel package in MATLAB, can significantly speed processing up  
    numcpu = 8;

% Field model (ONERA format) to load data for, and to then apply later for visualization
    %mfmtxt = 'T89'; % [T89], T96, T01, T01s, T04s, TS07Dmid15
    mfmtxt = 'TS07Dmid15'; % [T89], T96, T01, T01s, T04s, TS07Dmid15
        use_high_p = false; % use the high precision MFM computations? These may not exist for all spc yet
        use_lanlstar = false; %use lanlstar computed data instead of oneara lstar? can be quite 
                          %inaccurate, and relies on additional solar wind inputs, not only Kp.

% Data and Model errors for the Kalman Filter
    load_errors_from_file = false; % [true], false
        use_bias = false; % [false], true - utilize bias from analysis (not reliable atm)
        % true loads model errors and bias, and also data errors.
        %err_file_dir = 'reanalysis_files/';
%        err_file_dir = 'reanalysis_ts07d_err/';
        err_file_dir = 'reanalysis_t04s_err/';
        % If false then choose fixed errors to apply for model and data
        % if [true] then the errors will be loaded from a file specified in the Ini_dirs.m
        % Note that for best results with fixed errors, use equal errors [e.g. Kellerman et al., 2014]
        fixed_mod_err = 0.5; %fixed model error, percent of variance
        fixed_mod_bias = 0; %fixed model bias, as percent of PSD
        fixed_obs_err = 0.5; %fixed data error, percent of variance
%         fixed_mod_err = 0.1; %fixed model error, percent of variance
%         fixed_mod_bias = 0; %fixed model bias, as percent of PSD
%         fixed_obs_err = 500; %fixed data error, percent of variance      
        
        
% compute model forecast? it saves model psd alongside the 3DDA for later analysis
    save_verb_forecast = true; % true, [false] - save verb variable for error determination

% GREEP data may be loaded to assimilated if no spacecraft data are available.
    % uses Invariant GREEP or IGREEP at L* = 6
    assimilate_greep = false; % [false], true

%run VERB without assimilation? Still loads matrices each time step and computes solution in MATLAB
    model_simulation_only = false; % [false], true
        %Code doesn't work properly when false, need to fix variable scope
        run_cpp_only = false; % skip the matrices output/matlab solution and compute in cpp only
                             % this defaults to tridiag for the solution method although one could change 
                             % this at will

% The VERB code will be run from this directory, it will be created if it doesn't already exist
    %simulation_dir=['VERB2p3_',mfmtxt,'_err']; % directory, executables etc will be create/copied to the new folder
    simulation_dir=['VERB2p3_',mfmtxt,'_final2']; % directory, executables etc will be create/copied to the new folder
        % if you want to copy Dxx_ini and DiffCoeff files to the simulation dir then specify here
        copy_dxx_files = true;
            % Note that this only occurs if the folder doesn't already exist!
            source_simulation_dir = 'VERB/'; % specify the simulation base dir to copy from
                                                  % an error will be thrown if it doesn't exist

% Choose the solution method for inverting A to solve for f[t+1] in MATLAB
    % Note that our equation is "A*f[t+1] = B*f[t] + C"
    sol_meth = 'Gaussian'; %['Gaussian'], 'inverse','tridiag'
        % inverse - original solution, potentially inaccurate in matlab, see MATLAB doc's
        % tridiag - tridiagonal solution, Thomas algorithm
        % Gaussian - Guassian elimination (note that matrix A is tridiagonal anyway by default)

% capability to load diffusion coefficients and interpolate them to the grid we want
    use_non_field_aligned_Dxx = true; %[false], true
        % default [false] will compute the FA diffusion coefficients in C++
        load_Dxx_version = 3; %0 - FA, 1 - Ksenia Typical, 2 K14, 3 Dxx16, 4 BinBin

% simulate magnetopause shadowing.
    use_MP = false; % [true], false
        % simulated as LCDS from T89 with a centered dipole, and a lossy outer boundary in L in VERB

% set up grid structure
    %% define the VERB grid, need this for computation and reanalysis
    L_min = 1;
    L_max = 6.6;
    E_min = 0.01;
    E_max = 10;
    a_min = 0.3;
    a_max = 89.7;

    NL = 29;
    NE = 101;
    NA = 91;

% we can restart the simulation from a previous filedate, provided there is data (not implemented yet)
    use_SS = true; %[true], false 
        % default is to start from SS solution, however you can load a file if you want
        restart_file = 'PSD_Reanalysis_newest'; % name of file to start from
            % you can specify a time in the file to load the initial PSD from
            % default will load the last time in the file
            use_restart_time = false; %[false], true
                restart_time = datenum(2015,10,21,7,28,0);


% plotting settings if you want to visualize output as the assimilation runs
    plot_PSD_profiles = false; %plot PSD profiles across InvMu for a fixed alpha?
    prog_plot = false; % plot invariant flux for a chosen fixed pitch angle and energy?
        target_epc = 0.5;  %energy to plot
        target_alpha = 50; %pitch angle to plot

% Do you want to specify an output filename? 
    user_output_name = false; %true or [false],
        % The default [false] output is recommended, as it specifies a name representative
        % of the options you choose below, and the time period
        output_file_name = 'some_name.mat'; %if above switch is true then specify your name here

% specify the output file directory
    %output_dir = 'reanalysis_ts07d_err/'; % will be created if it doesn't exist
        % note that this will be overwritten if you make multiple simulations in the run script
    
    output_dir= ['reanalysis_',mfmtxt,'_final'];

    use_matr_taulc = false; % [false], true. Uses the version of the code with losses built into matr_A
                            % alternative is to use taulc from sources and losses.
                            % the latter matches paper results

    use_latest_VERB = true;     % [true] ,false - if true then the above switch will not work
                            % uses latest VERB 2.# version

    addpath([getenv('HOME'),'/Documents/VERB/22. 3D VERB DA_Spectrum_LT/Code/Various_functions/']);

    % load file with plasmapause location
    Lpp = load('Parameters/Lpp_01Oct15_01Oct16.txt');
    
    % load file with dynamic pressure
    Pdyn = load('Parameters/pdyn_01Oct15_01Oct16.txt');

    % load file with Bf
%     Bf = load('Parameters/BF_01Mar13_31Mar13_interp.dat');
    
    % load file with magnetopause location
    Mp = load('Parameters/mp_01Oct15_01Oct16.txt');
    
% load file with Kp
%     Kp = load('Parameters/kp_01Mar13_31Mar13.dat');
    dummy = zeros(1,92);    
    
    it_tot = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% end settings 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
