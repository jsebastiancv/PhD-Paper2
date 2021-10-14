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

%% This section is dedicated to setting up the output file names, loading appropriate data etc.

output_dir = [output_dir,'/']; %mark it
output_dir = regexprep(output_dir,'//','/'); %fix double mark

d_time = timestep/24 ; %time step in days

if use_lanlstar
    fname='lanl';
else
    fname='onera';
end

if model_simulation_only
    fname='No_DA';
else
    run_cpp_only = false; %default 
end

% loads magnetopause data, only T89 ATM. We could implement a computation for the 
% period of interest, though I think this would be computationally expensive for the time being
if use_MP
    kpdata = load(['files/MP_centeredDip__LCDS_T89.txt']);
    mp_data = reshape(kpdata,18,27,4);
    mp_Kp = squeeze(mp_data(1,:,1));
end

if ~run_cpp_only % no need if running cpp
    try
        parpool(numcpu);
    catch end
end

time_range = [start_time end_time];
if ~user_output_name
    if model_simulation_only
        if run_cpp_only
            sftxt = 'CPP_VERBonly';
        else
            sftxt = 'MS_VERBonly';
        end
        if use_MP
            sftxt = [sftxt,'_MP'];
        else
            sftxt = [sftxt,'_noMP'];
        end
        sftxt = [sftxt,'_',yrmo];
        if run_cpp_only
            output_file_name = ['./Output/',sftxt,'_',fname,'_',mfmtxt,'.mat'];
        else
            output_file_name = ['./Output/',sftxt,'_',sol_meth,'_',fname,'_',mfmtxt,'.mat'];
        end
    else
        % output file string creation
        sftxt = 'Reanalysis';
        if assimilate_greep
            sftxt = [sftxt,'_GREEP'];
        end
        if use_MP
            sftxt = [sftxt,'_MP'];
        else
            sftxt = [sftxt,'_noMP'];
        end
        sftxt = [sftxt,'_',yrmo];
        output_file_name = [output_dir,sftxt,'_',sol_meth,'_',fname,'_',mfmtxt,'.mat'];
    end
    VERB_suff = '';

    if use_latest_VERB
        output_file_name = regexprep(output_file_name,'Reanalysis_','Reanalysis_LatestVERB_');        
        VERB_exe = 'VERB_code_2.3';
    else
        VERB_exe = 'VERB_code_2.0';
        if ~use_matr_taulc
            output_file_name = regexprep(output_file_name,'Reanalysis_','Reanalysis_SLTauLC_');
            VERB_suff = '_ImpLoss';
        end
    end
    if ~load_errors_from_file
        output_file_name = regexprep(output_file_name,'Reanalysis_','Reanalysis_EQE_');
    end
    if save_verb_forecast
        output_errfile_name = regexprep(output_file_name,'Reanalysis','ReanalysisErr');
    end
    [~,name,ext] = fileparts(output_file_name);
    file_string = [name,ext];

    [~,~] = mkdir(output_dir,'f');
    % note that currently the mfmtxt doesn't mean anything for VERB only simulations
    % in future it will point to an MFM dependent MP in L.
end

%% here we check of the output file already exists, and load it if necessary. This will skip the loop.
fck = dir(output_file_name);
if ~isempty(fck)
    fprintf('reanalysis file found, loading...');
    use_SS = false; %false after the first instance
    load(output_file_name);
    last_PSD = squeeze(SimPSD(end,:,:,:));
    fprintf('done\n');
else
    % just for plotting
    target_pc = pfunc(target_epc);
    target_alpha = degtorad(target_alpha);

    SimTime = time_range(1):d_time:time_range(end);
    nt = length(SimTime);
    hist_kp

    % set up Kp, remove any NaN's
    % since Kp represents activity over time, we should interpolate to the time point
    % midway between each time step: if simulating from t = 0 to t = 1 then we should 
    % find the kp at t = 0.5 to best represent the interval.

    % the variable simulation time represents the assimilation of data on time step t = 1
    % and model results from time step t = 0 to t = 1; 
    % therefore, Kp should be for t = t-dt/2;
    % since we run the simulation one index behind (start at index t, but reference index t-1), we should actually
    % compute kp for t = t+dt/2 to represent the simulation from t=1 to t=2;

    Kp_t = SimTime + d_time/2;
    zz=find(~isnan(kp_index));
    Kp = interp1(kp_time(zz),kp_index(zz),Kp_t,'nearest','extrap');
    Kp_recent = -1;  % we introduce this factor for the case when we load precomputed ABC matrices (suboptimal matrices, but faster. That method is not particularly accurate but may give a first approximation)
   
    if use_SS %use steady state solution as IPSD
        Kp_sim = Kp(1);
        VERB_SS %run simulation to get initial PSD
        folder = [simulation_dir,'Output/'];
        [L, epc, alpha, pc] = load_plt([folder, 'perp_grid.plt']);
        SimL = L.arr; SimInvEnergy = epc.arr; SimInvAlpha = alpha.arr; pc = pc.arr;

        %specification of initial analysis error covariance matrices
        if ~exist('Pa_alpha','var');
            Pa_alpha=zeros(NL,NE,NA,NA); % then initialize as zero
            Pa_pc=zeros(NL,NA,NE,NE);
            Pa_L=zeros(NE,NA,NL,NL);
        end

        SimInvMu = pc2mu(SimL,pc,SimInvAlpha);
        SimInvK = Lalpha2K(SimL,SimInvAlpha);
        use_SS = false; %false after the first instance
    else
        IPSD = last_PSD; %
    end
    IPSD(IPSD <= 1e-21) = 1e-21;

    folder = [simulation_dir,'Output/'];
    [L, epc, alpha, pc] = load_plt([folder,'perp_grid.plt']);
    SimL = L.arr; SimInvEnergy = epc.arr; SimInvAlpha = alpha.arr; pc = pc.arr;
    SimInvMu = pc2mu(SimL,pc,SimInvAlpha);
    SimInvK = Lalpha2K(SimL,SimInvAlpha);
    ModPSD(1,:,:,:) = IPSD;
    
    L_vec = squeeze(L.arr(:,1,1));
    
    l_idx = find(SimL(:,1,1) == 6.6); % index for L* = 6.6
    [~,e_idx] = min(abs(SimInvEnergy(l_idx,:,1) - 0.8 )); % index for E closest to 1 MeV at L* = 6.6;

    angle = 90 * pi/180;
    [~,a_idx] = min(abs(SimInvAlpha(l_idx,e_idx,:) - angle)); % index for alpha closest to 90 degree at E = 1 MeV and L* = 6.6;

    PSD_init = ModPSD(1,l_idx,e_idx,a_idx); % we will use this value to rescale Bf in the one-step code

    lvals = SimL(:,end,end); %copied from line 138
       
    if ini_model_err
        model_error = ones(NL,NE,NA) .* fixed_mod_err; %default model error
        model_bias = ones(NL,NE,NA) .* fixed_mod_bias;
        ini_model_err = false;
    else
        if load_errors_from_file
            err_file = [err_file_dir,mfmtxt,'_errors.mat']; %VERB error file           
            fck = dir(err_file);
            if ~isempty(fck)
                load(err_file)
            else
                error('error file does not exist');
            end
            fprintf('gridding errors to current grid and time step...');
            m_err = griddata(reshape(errL,[],1),log10(reshape(errInvEnergy,[],1)),...
                reshape(errInvAlpha,[],1),reshape(verb_std,[],1),reshape(SimL,[],1),...
                reshape(log10(SimInvEnergy),[],1),reshape(SimInvAlpha,[],1),'natural');
            m_err(isnan(m_err)) = 0.5;
            m_err = m_err .* timestep; %errors are saved in hourly format, so scaling is trivial
            model_error = reshape(m_err,NL,NE,NA);
            fprintf('done\n');
            if use_bias
                fprintf('gridding bias to current grid and time step...');
                m_bias = griddata(reshape(errL,[],1),log10(reshape(errInvEnergy,[],1)),...
                    reshape(errInvAlpha,[],1),reshape(verb_bias,[],1),reshape(SimL,[],1),...
                    reshape(log10(SimInvEnergy),[],1),reshape(SimInvAlpha,[],1),'natural');
                m_bias(isnan(m_err)) = 1.0; %default 1.0, which is no bias
                m_bias = m_bias .* timestep; %bias is also saved in hourly format, so scaling is trivial
                model_bias = reshape(m_bias,NL,NE,NA);
                fprintf('done\n');
            else
                model_bias = ones(NL,NE,NA);
            end
            clear verb_std verb_bias m_err m_bias
        else
            model_error = ones(NL,NE,NA) .* fixed_mod_err;
            model_bias = ones(NL,NE,NA) .* fixed_mod_bias;
        end
    end

    if use_non_field_aligned_Dxx
        % then we must load the diffusion coefficients, interpolate them, then save them
        % in the new grid size
        Dxx_load_folder = [Satellite_data_dir,'DiffusionCoefficient/Matrices/'];
        Dxx_output_folder1 = [simulation_dir,'Dxx_ini/'];
        Dxx_output_folder2 = [simulation_dir,'DiffCoeff/'];
%         delete([Dxx_output_folder1,'*']);
%         delete([Dxx_output_folder2,'*']);

        if load_Dxx_version == 1 %Binbin
            load_folder = [Dxx_load_folder,'Binbin_Version/VERB2/'];
        elseif load_Dxx_version == 2 % Ksenia 14
            load_folder = [Dxx_load_folder,'Ksenia_Version/VERB2/'];
        end
    end

    %% we only need to load data if we are assimilating it
    if ~model_simulation_only
        if use_high_p
            hptxt = 'HPREC';
        else
            hptxt = ''; 
        end
        rb_data_dir = ['/export/.automounter/mag5/data/rbm-data/RBSP/']; 
        goes_data_dir = ['/export/.automounter/mag5/data/rbm-data/GOES/']; 
        
        sd = time_range(1);
        ed = time_range(2);

        spcs = {'rbspa','rbspb','rbspa','rbspb','GOES13','GOES15'};
        instrs = {'mageis','mageis','rept','rept','MAGEDandEPEAD','MAGEDandEPEAD'};
        spc_dirs = {rb_data_dir,rb_data_dir,rb_data_dir,rb_data_dir,goes_data_dir,goes_data_dir};
        
        nspc = length(spcs);
        if assimilate_greep
            % then we don't assimilate data...we can change this later
            nspc = 0;
            PSD_obs = nan(1,nt,NL,NE,NA);
        else
            PSD_obs = nan(nspc,nt,NL,NE,NA);
        end
        PSD_obs_err = nan(size(PSD_obs,1),nt,NL,NE,NA);
        for ispc=1:nspc
            cdate = sd;
            spc =spcs{ispc};
            sdir = [spc_dirs{ispc},spc,'/Processed_Mat_Files/'];
            instr = instrs{ispc};
            fprintf('loading %s %s...',spc,instr);
            try
                spc_time = [];
                spc_PSD = [];
                spc_errPSD = [];
                spc_Lstar = [];
                spc_InvK = [];
                spc_InvMu = [];
                while cdate < ed
                    [yr,mo,~,~,~,~] = datevec(cdate);
                    dates = make_dates(yr,mo);
                    if strcmp(instr,'MAGEDandEPEAD'); file = [sdir,spc,'_',instr,'_',dates,'_psd_n4_4_',hptxt,mfmtxt,'_ver4_calibrated.mat']; end
                    if strcmp(instr,'mageis'); file = [sdir,spc,'_',instr,'_',dates,'_psd_n4_4_',hptxt,mfmtxt,'_ver4_calibrated.mat']; end
                    if strcmp(instr,'rept'); file = [sdir,spc,'_',instr,'_',dates,'_psd_',hptxt,'ver4.mat']; end
                    if load_errors_from_file
                        try
                            load(file);
                            file = [sdir,spc,'_',instr,'_',dates,'_errpsd_',hptxt,...
                                mfmtxt,'_ver4.mat'];
                        catch
                            error('No error file for %s on %s\n',spc, cdate);
                        end
                    end
                    load(file);
                    
                    if strcmp(instr,'MAGEDandEPEAD'); file = [sdir,spc,'_',instr,'_',dates,'_lstar_',hptxt,mfmtxt,'_ver4.mat']; end
                    if strcmp(instr,'mageis'); file = [sdir,spc,'_',instr,'_',dates,'_lstar_',hptxt,mfmtxt,'_ver4.mat']; end
                    if strcmp(instr,'rept'); file = [sdir,spc,'_',instr,'_',dates,'_lstar_',hptxt,mfmtxt,'_ver4.mat']; end                    
                    load(file);
                    
                    if strcmp(instr,'MAGEDandEPEAD'); file = [sdir,spc,'_',instr,'_',dates,'_invmu_and_invk_',hptxt,mfmtxt,'_ver4.mat']; end
                    if strcmp(instr,'mageis'); file = [sdir,spc,'_',instr,'_',dates,'_invmu_and_invk_',hptxt,mfmtxt,'_ver4.mat']; end
                    if strcmp(instr,'rept'); file = [sdir,spc,'_',instr,'_',dates,'_invmu_and_invk_',hptxt,mfmtxt,'_ver4.mat']; end                           
                    load(file);

                    if strcmp(instr,'MAGEDandEPEAD')
                        file = [sdir,spc,'_',instr,'_',dates,'_alpha_and_energy_',hptxt,mfmtxt,'_ver4.mat']; 
                        load(file)
                        ok_ech=find(energy_channels <= 2); 
                        PSD = PSD(:,ok_ech,:);
                        InvMu = InvMu(:,ok_ech,:);
        %                errPSD = errPSD(:,ok_ech,:);
                    elseif strcmp(instr,'mageis')
                        file = [sdir,spc,'_',instr,'_',dates,'_alpha_and_energy_n4_4_',hptxt,mfmtxt,'_ver4.mat']; load(file);
                        for itm=1:length(time)
                            bad_ech=find(energy_channels(itm,:) > 2);  
                            if ~isempty(bad_ech)
                                PSD(itm,bad_ech,:) = NaN;    
                                InvMu(itm,bad_ech,:) = NaN;  
                                errPSD(itm,bad_ech,:) = NaN;
                            end 
                        end 
                    end    
  

                    spc_time = cat(1,spc_time,time);
                    spc_PSD = cat(1,spc_PSD,PSD);
                    if load_errors_from_file
                        spc_errPSD = cat(1,spc_errPSD,errPSD);
                    end
                    spc_Lstar = cat(1,spc_Lstar,Lstar);
                    spc_InvK = cat(1,spc_InvK,InvK);
                    spc_InvMu = cat(1,spc_InvMu,InvMu);
                    dim = find_DIM(yr,mo);
                    cdate = datenum(yr,mo,dim) + 1;
                end
                fprintf('interpolating %s %s PSD:\n',spc,instr);
                [PSD_g] = interp2grid(SimTime,SimInvMu,SimInvK,SimL,...
                    spc_time,spc_InvMu,spc_InvK,spc_Lstar,spc_PSD);    
                PSD_obs(ispc,:,:,:,:) = PSD_g .* 2.997e7 ; %convert to VERB units
                clear PSD_g

                if load_errors_from_file
                    fprintf('interpolating %s %s PSDerr:\n',spc,instr);
                    [PSDerr_g] = interp2grid(SimTime,SimInvMu,SimInvK,SimL,...
                        spc_time,spc_InvMu,spc_InvK,spc_Lstar,spc_errPSD);    
                        %errPSD is a multiplication factor
                    PSD_obs_err(ispc,:,:,:,:) = PSDerr_g;
                    clear spc_errPSD PSD_err_g
                
                end
                clear spc_time spc_PSD spc_Lstar spc_InvK spc_InvMu

            catch
                fprintf('no data for %s %s\n',spc,instr)
            end
        end
                
        if (assimilate_greep)
            fprintf('loading GREEP data...');
            greep_mfm = 'T89';

            greepdir = ['/home/akellerman/projects/GREEP_code/branches/invariant_greep/matlab/mat/'];
            cdate = start_time;
            greep_time = [];
            greep_PSD = [];
            greep_InvK = [];
            greep_InvMu = [];
            greep_Lstar = 6;

            while cdate < end_time
                [yr,~,~,~,~,~] = datevec(cdate);
                greep_file = [greepdir,'greep_hist_forecast_',num2str(yr),'_',greep_mfm,'.mat'];
                load(greep_file); 
                greep_InvMu = muvals(2:end-1);
                greep_InvK = kvals;
                greep_PSD = cat(1,greep_PSD,forecast_PSD(:,2:end-1,:));
                greep_time = cat(1,greep_time,forecast_time);
                clear muvals kvals forecast_PSD forecast_time
                cdate = datenum(yr,12,31)+1;
            end
            zz=find(greep_time >= sd & greep_time <= ed);
            if ~isempty(zz)
                greep_PSD = greep_PSD(zz,:,:);
                greep_time = greep_time(zz);
                greep_Lstar = 6;
                fprintf('done\n');
                fprintf('interpolating GREEP forecast:\n');
                [greep_PSD_g] = interpGREEP2grid(SimTime,SimInvMu,SimInvK,SimL,greep_time,greep_InvMu,...
                    greep_InvK,greep_Lstar,greep_PSD);
                clear greep_PSD;
                PSD_obs(end,:,:,:,:) = greep_PSD_g .* 2.997e7 ; %convert to VERB units
                if load_errors_from_file
                    PSD_obs_err(end,:,:,:,:) = greep_PSD_g.*0 + 0.5; % to do, actually determine the errors, it can be done
                    clear greep_PSD_g
                end
                clear greep_time greep_InvMu greep_InvK greep_Lstar greep_PSD
            else
                assimilate_greep = false;
                PSD_obs = PSD_obs(1:end-1,:,:,:,:);
                if load_errors_from_file
                    PSD_obs_err = PSD_obs_err(1:end-1,:,:,:,:);
                end
            end
        end
    end

    %% obs errors
    % errors = sum ( w.*PSD) ./ (w)
   if load_errors_from_file
        fprintf('combining errors and PSD...');
        t1 = tic;
        zz = find(isnan(PSD_obs));
        PSD_obs_err(zz) = NaN;
        num = nansum(PSD_obs ./ PSD_obs_err,1);
        den = nansum(1./PSD_obs_err,1) ; % PSD_obs may have nans, while w may not
        PSD_obs = squeeze(num ./ den);
        clear num den
        PSD_obs_err = squeeze(nanmean(PSD_obs_err,1));
        errt = toc(t1);
        [tstr] = format_time(errt);
        fprintf('done in %s\n',tstr); 

        % removing nans from bias
        for ie=1:NE
            for ia=1:NA
                zz=find(~isnan(model_bias(:,ie,ia)));
                if ~isempty(zz)
                    model_bias(:,ie,ia) = interp1(SimL(zz,ie,ia),model_bias(zz,ie,ia),...
                        SimL(:,ie,ia),'nearest','extrap');
                end
            end
        end
    else
        PSD_obs = squeeze(nanmean(PSD_obs,1));
        PSD_obs_err = ones(NL,NE,NA) .* fixed_obs_err;        
    end

    SimPSD = nan(nt,NL,NE,NA);
    if save_verb_forecast
        ModPSD = nan(nt,NL,NE,NA);
        ModPSDup = zeros(nt,NL,NE,NA);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %%%     Start simulation/assimilation 
    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

    %% you can run the C++ code simulation only, this skips the for-loop below all together
    if run_cpp_only 
        % then we can initialize the computation in cpp only and skip the for loop below.

        % In this case, we want to interpolate Kp to the mid time of the simulation interval, 
        % we do this, as VERB will read the kp and update Kp
        % by interpolating to the closest kp points in time surrounding the current simulation time point
        % the current simulation time represents the END of the interval, so really we need to shift the Kp by 
        % dt/2. 

        % note that the simulation time here is similar to that for the data assimilation, but the VERB code
        % referencing is different. VERB references index t for time interval t so we must set the Kp time 
        % differently

        % Simulation interval
        % t=0 -------dt-------t=1-------dt-------t=2;

        % For simulation at t=1 (representative of 0 to 1) VERB chooses Kp value at Kp time interpolated linearly
        % to t=1;

        % Hence to make sure Kp is accurate, we should interpolate Kp appropriatly

        % in Load_kp we load kp_time as the midpoint in the kp interval (1.5 represents 0-3 hrs)
        % so to interpolate it to t=1 we should interpolate the kp value to (t=0+dt/2)
        % this is the same as what we have done above for reanalysis, so we can keep the Kp

        % Then to read it correctly in VERB we should adjust the time to equal the SimTime + dt;
        % this way, the kp read into VERB is representative of the interval which is simulated 
        % and not interpolated in some strange way. This also forces a ~= 0 in the Parameters.cpp, load_1d;

        Kp_tverb = SimTime + d_time - SimTime(1); %sim starts from 0 

        if size(Kp_tverb,1) == 1;
            Kp_tverb = Kp_tverb';
        end
        if size(Kp,1) == 1;
            Kp = Kp';
        end
        timekp = cat(2,Kp_tverb,Kp);
        saveascii(timekp,[simulation_dir,'Input/kp.dat'],5)
        ndays = time_range(2) - time_range(1);

        output_step = 1; %each step = 1

        mp_set = false; %no mp yet

        lastPSD.time = 0;
        lastPSD.arr = squeeze(IPSD);  %only if we load from a file
        
        VERB_noMM_new 
        %new verb_noMM, just calls new executable with solveMatrix
        %parameter = false
        
        SimPSD = PSD.arr;
        save(output_file_name,'SimTime',...
            'SimPSD','Kp');

        return
    end

    SimPSD(1,:,:,:) = IPSD; % set the first point to IPSD
    %ModPSD(1,:,:,:) = NaN; % set the first point to NaN, it is already so no need to do this

    lastPSD.time = 0;
    lastPSD.arr = IPSD; % PSD that is given to VERB for the first time step

    fprintf('\n');
    fprintf('Beginning Simulation\n\n');
    totaltime = tic;
    
%     L_up_initial(:,:) = SimPSD(1,end,:,:); % added
    pc_low_initial(:,:) = SimPSD(1,:,1,:); % added
    
%     it_month = it_month + 1;
    
    for it = 2 : nt % simulation always starts at point 2 (first computation is from f(t[1]) to f(t[2]);
        
        it_tot = it_tot + 1;
        
        fprintf('### 3D Data Assimilative VERB code ###\n')
        fprintf('Out file: %s\n',file_string)
        loopt = tic;
        Kp_sim = Kp(it); % what was the kp representative of the transition of PSD
                            % to the current time? We will use this in the simulation
        
        %% for MP stuff, find the corresponding Kp value
%         if use_MP
%             kpi=find(mp_Kp == round(Kp_sim*10));
%             mp_L = squeeze(mp_data(:,kpi,3));
%             mp_K = squeeze(mp_data(:,kpi,4));
%         end
%         minL=min(mp_L);
        
        fprintf('Time: %s\n', datestr(SimTime(it),'yyyy-mm-dd, HH:MM'));
        fprintf('Step: %s/%s Kp: %s\n',num2str(it),num2str(nt),num2str(Kp_sim,3));

        if mod(it,10) == 0 && plot_PSD_profiles
            yrng = [-21 10];
            trng = [1 7];
            subplot(4,2,1)
            plot(SimL(:,1,1),squeeze(log10(SimPSD(it-1,:,:,12))));
            xlim(trng);
            ylim(yrng);

            subplot(4,2,2)
            plot(SimL(:,1,1),squeeze(log10(PSD_obs(it,:,:,12))));
            xlim(trng);
            ylim(yrng);
        end

        %% introduce loop-specific variables
        % The observations are representative of the final state,
        % or an average over the past several steps if desired, f[t+1];
        
        if it < 5
            PSD_d_loop = squeeze(nanmean(PSD_obs(1:it,:,:,:),1));
        else
            PSD_d_loop = squeeze(nanmean(PSD_obs(it-(ndastep-1):it,:,:,:),1));
        end
      
        if load_errors_from_file
            PSD_derr_loop = squeeze(nanmean(PSD_obs_err(it-(ndastep-1):it,:,:,:),1));
        else
            PSD_derr_loop = PSD_obs_err;
        end
        PSD_d_loop(PSD_d_loop <= 1e-21) = 1e-21;

        % reanalysis/simulation
        % the simulation begins at f[t] and we compute f[t+1], so inx is t-1 here
        PSD_r_loop = squeeze(SimPSD(it-1,:,:,:));
        PSD_r_loop(PSD_r_loop <= 1e-21 | isnan(PSD_r_loop)) = 1e-21;



        %% run VERB for this time step
        lastPSD.arr = PSD_r_loop;
%         lastPSD.arr(end,:,:) = L_up_initial; % added
        lastPSD.arr(:,1,:) = pc_low_initial; % added
           
%         Bf_sim = Bf(it,2);
%         
%         Bf_sim_rec = Bf_sim * (PSD_init / (lastPSD.arr(l_idx,e_idx,a_idx)));        
%         
        Lpp_sim = Lpp(it_tot,2);
        
        Pdyn_sim = Pdyn(it_tot,2);        
%         
%         if Pdyn_sim > 50
%             Pdyn_sim = 0;
%         end
        
        save_spectrum

        Mp_it = Mp(it_tot,:);
        Mp_it(1) = 0;% 0.0416667;

        extra = Mp_it;
        extra(1) = 0.0416667;

        min_mp = min(Mp_it(2:end));
        
        %% 1. Forecast performed with one-step VERB (3D + Mixed terms) 

        if use_MP
            minL = min(mp_L);
            if minL < lvals(end) % then load MP coefficients for time step
                mp_set = true;
            else
                mp_set = false;
            end            
        else
            mp_set = false;
        end    
        
        VERB_one_step
        
        dataFileName = [simulation_dir,'/Output/OutPSD.dat'];   
        skip_zones = 0;
        PSD_MT = load_plt(dataFileName, 'permute', 'n_zones', 3000, 'skip_zones', skip_zones);      
        PSD_m_loop = squeeze(PSD_MT.arr(2,:,:,:)); %% PSD_m_loop = one-step prediction   
        PSD_m_loop(PSD_m_loop <= 1e-21 | isnan(PSD_m_loop)) = 1e-21;
%         
%         for ia = 1:91
%             for il = 1:29
%                 if  L_vec(il) > Mp_it(1,ia+1) 
%                     PSD_m_loop(il,:,ia) = 1e-21;
%                 end
%             end
%         end        
                
        %% 2. Assimilation / update using split-operator technique    
%         
        if Pdyn_sim > 3 
             if min_mp < 6.6
                 mat_file = ['Model Matrices/Debug_MT_Val_EMICs_waverate4%/Kp=',num2str(Kp_sim),'.mat'];
             else
                mat_file = ['Model Matrices/Debug_MT_Der_EMICs_waverate4%/Kp=',num2str(Kp_sim),'.mat'];
             end
        else
             if min_mp < 6.6
                 mat_file = ['Model Matrices/Debug_MT_Val_NoEMICs/Kp=',num2str(Kp_sim),'.mat'];
             else
                mat_file = ['Model Matrices/Debug_MT_Der_NoEMICs/Kp=',num2str(Kp_sim),'.mat'];
             end
        end
%             mat_file = ['Model Matrices/DiffMat_Hiss2016_Chorus2017_NoEMICs_ZeroDer/Kp=',num2str(Kp_sim),'.mat'];
            load(mat_file);


        assimt = tic;        
        %%diffusion for L
%         KD_L = zeros(29,101,91);
        parfor j=1:NE
            PSD_L=PSD_m_loop(:,j,:);
%             Kd_L = zeros(size(PSD_L));%(29,1,91);
            P_L = Pa_L(j,:,:,:);
            for k=1:NA
                % get model operator
                [~,matmodel]=matsolve(L_matr_A(j,k,:,:),...
                    L_matr_B(j,k,:),L_matr_C(j,k,:),PSD_L(:,1,k),'solution_method',sol_meth);

                if (model_simulation_only==0 & mod(it,ndastep)==0)
                    % assimilate data
                  [PSD_L(:,1,k),P_L(1,k,:,:)] = KalmanFilter_1Derrs(...
                        squeeze(PSD_L(:,1,k)),squeeze(model_error(:,j,k)),matmodel,...
                        squeeze(P_L(1,k,:,:)),squeeze(PSD_d_loop(:,j,k))',...
                        squeeze(PSD_derr_loop(:,j,k)));
                end
            end
            PSD_r_loop(:,j,:) = PSD_L;
            Pa_L(j,:,:,:) = P_L;
%             KD_L(:,j,:) = Kd_L;
        end
        
%         if mp_set % comment this is MP = false
%             if minL < L_max
%                 for iL=1:NL
%                     if lvals(iL) >= minL
%                         for iK = 1:NA
%                             Kval = SimInvK(iL,iK);
%                             max_L = interp1(mp_K,mp_L,Kval,'nearest','extrap');
%                             if lvals(iL) >= max_L %then we are outside 
%                                 PSD_r_loop(iL,:,iK) = 1e-21; % set to zero
% %                                 if save_verb_forecast
% %                                     PSD_m_loop(iL,:,iK) = 1e-21; % set to zero
% %                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
        

%         for ia = 1:91
%             for il = 1:29
%                 if  L_vec(il) > Mp_it(1,ia+1) 
%                     PSD_r_loop(il,:,ia) = 1e-21;
%                 end
%             end
%         end        
                
        PSD_r_loop(PSD_r_loop <= 1e-21 | isnan(PSD_r_loop)) = 1e-21;
        
        %%diffusion for alpha
%         KD_alpha = zeros(29,101,91); 
        parfor j=1:NE
            PSD_A=PSD_r_loop(:,j,:); 
%             Kd_alpha = zeros(size(PSD_A));
            P_al = Pa_alpha(:,j,:,:);
            for i=1:NL
                % get model operator
                [~,matmodel]=matsolve(alpha_matr_A(i,j,:,:),...
                    alpha_matr_B(i,j,:),alpha_matr_C(i,j,:),PSD_A(i,1,:),'solution_method',sol_meth);  

                if (model_simulation_only==0 & mod(it,ndastep)==0)
                    % assimilate data
                    [PSD_A(i,1,:),P_al(i,1,:,:)] = KalmanFilter_1Derrs(...
                        squeeze(PSD_A(i,1,:)),squeeze(model_error(i,j,:)),matmodel,...
                        squeeze(P_al(i,1,:,:)),squeeze(PSD_d_loop(i,j,:))',...
                        squeeze(PSD_derr_loop(i,j,:)));
                end
            end
            PSD_r_loop(:,j,:) = PSD_A;
            Pa_alpha(:,j,:,:) = P_al;
%             KD_alpha(:,j,:) = Kd_alpha;
        end

% %         if mp_set % comment this is MP = false
% %             if minL < L_max
% %                 for iL=1:NL
% %                     if lvals(iL) >= minL
% %                         for iK = 1:NA
% %                             Kval = SimInvK(iL,iK);
% %                             max_L = interp1(mp_K,mp_L,Kval,'nearest','extrap');
% %                             if lvals(iL) >= max_L %then we are outside 
% %                                 PSD_r_loop(iL,:,iK) = 1e-21; % set to zero
% % %                                 if save_verb_forecast
% % %                                     PSD_m_loop(iL,:,iK) = 1e-21; % set to zero
% % %                                 end
% %                             end
% %                         end
% %                     end
% %                 end
% %             end
% %         end

%         for ia = 1:91
%             for il = 1:29
%                 if  L_vec(il) > Mp_it(1,ia+1) 
%                     PSD_r_loop(il,:,ia) = 1e-21;
%                 end
%             end
%         end        
              
        PSD_r_loop(PSD_r_loop <= 1e-21 | isnan(PSD_r_loop)) = 1e-21;

        %%diffusion for energy
%         KD_pc = zeros(29,101,91);
        parfor k=1:NA
            PSD_E = PSD_r_loop(:,:,k); 
%             Kd_pc = zeros(size(PSD_E));
            P_pc = Pa_pc(:,k,:,:);
            for i=1:NL
                % get model operator
                [~,matmodel]=matsolve(pc_matr_A(i,k,:,:),...
                    pc_matr_B(i,k,:),pc_matr_C(i,k,:),PSD_E(i,:,1),'solution_method',sol_meth);

                if (model_simulation_only==0 & mod(it,ndastep)==0) % assimilate data
                    [PSD_E(i,:,1),P_pc(i,1,:,:)] = KalmanFilter_1Derrs(...
                        squeeze(PSD_E(i,:,1)),squeeze(model_error(i,:,k))',matmodel,...
                        squeeze(P_pc(i,1,:,:)),squeeze(PSD_d_loop(i,:,k))',...
                        squeeze(PSD_derr_loop(i,:,k))');
                end 
            end                           
            PSD_r_loop(:,:,k) = PSD_E;
            Pa_pc(:,k,:,:) = P_pc;
%             KD_pc(:,:,k) = Kd_pc;
        end
%  
% %         if mp_set % comment this is MP = false
% %             if minL < L_max
% %                 for iL=1:NL
% %                     if lvals(iL) >= minL
% %                         for iK = 1:NA
% %                             Kval = SimInvK(iL,iK);
% %                             max_L = interp1(mp_K,mp_L,Kval,'nearest','extrap');
% %                             if lvals(iL) >= max_L %then we are outside 
% %                                 PSD_r_loop(iL,:,iK) = 1e-21; % set to zero
% % %                                 if save_verb_forecast
% % %                                     PSD_m_loop(iL,:,iK) = 1e-21; % set to zero
% % %                                 end
% %                             end
% %                         end
% %                     end
% %                 end
% %             end
% %         end
% % 
%         for ia = 1:91
%             for il = 1:29
%                 if  L_vec(il) > Mp_it(1,ia+1) 
%                     PSD_r_loop(il,:,ia) = 1e-21;
%                 end
%             end
%         end        
               
         PSD_r_loop(PSD_r_loop <= 1e-21 | isnan(PSD_r_loop)) = 1e-21;

 %% ## END assimilation and begin plotting/populating arrays ##
        
          % save loop variables back to final arrays
        SimPSD(it,:,:,:) = PSD_r_loop; % Xa and Xf mix --- output from Kalman filter
%         KD_L_final(it,:,:,:) = KD_L;
%         KD_pc_final(it,:,:,:) = KD_pc;
%         KD_alpha_final(it,:,:,:) = KD_alpha;
                
        if save_verb_forecast
           updatepts = find(~isnan(PSD_d_loop));
           ModPSD(it,:,:,:) = PSD_m_loop; % pure Xf --- output from one-step VERB
           PSDup = zeros(NL,NE,NA);
           PSDup(updatepts) = 1;
           ModPSDup(it,:,:,:) = PSDup;
        end

        % some summary/progress information including an estimated completion time
        tstr = format_time(toc(assimt));
        fprintf('ASSIMILATION TIME: %s\n',tstr);
        tstr = format_time(toc(loopt));
        fprintf('TOTAL STEP TIME: %s\n',tstr);
        tstr = format_time(toc(totaltime));
        fprintf('TOTAL SIMULATION TIME: %s\n',tstr);
%         fprintf('Estimated Completion: %s\n',estimation_time(totaltime,it-1,nt-1));
        fprintf(' \n')
    end

    %% we are finished with the time period, now save our arrays

    % this is used for restarting the simulation if we are computing several months
    last_PSD = squeeze(SimPSD(end,:,:,:));

    % save output
    save(output_file_name,'SimPSD','SimTime','Kp','SimInvMu','SimInvK','SimInvEnergy','SimInvAlpha','SimL','Pa_alpha','Pa_pc','Pa_L','PSD_obs');

    % if computing Xf then we save it here
    if save_verb_forecast
        save(output_errfile_name,'ModPSD','ModPSDup');
    end
end
