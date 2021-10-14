%% parameters structure:
    steps=1;
    fid = fopen([simulation_dir,'Parameters.ini'], 'w');
    fprintf(fid, ['timeStep = ', num2str(timestep), '\n']);
    fprintf(fid, ['nDays = ', num2str(d_time), '\n']);
    fprintf(fid, ['general_Output_parameters.iterStep = ',num2str(steps), '\n']);
    
    fprintf(fid, ['max_threads = 20, \n']);
    
    fprintf(fid, 'useKp = constant\n');
    fprintf(fid, ['constKp = 2\n']);
        
    fprintf(fid, 'useBf =  No\n');
    
    fprintf(fid, ['localDiffusionsGrid_L.size = ', num2str(NL), '\n']);
    fprintf(fid, ['localDiffusionsGrid_L.min = ', num2str(L_min), '\n']);
    fprintf(fid, ['localDiffusionsGrid_L.max = ', num2str(L_max), '\n']);

    fprintf(fid, ['localDiffusionsGrid_epc.size = ', num2str(NE), '\n']);
    fprintf(fid, ['localDiffusionsGrid_epc.min = ', num2str(E_min), '\n']);
    fprintf(fid, ['localDiffusionsGrid_epc.max = ', num2str(E_max), '\n']);

    fprintf(fid, ['localDiffusionsGrid_alpha.size = ', num2str(NA), '\n']);
    fprintf(fid, ['localDiffusionsGrid_alpha.min = ', num2str(a_min), '\n']);
    fprintf(fid, ['localDiffusionsGrid_alpha.max = ', num2str(a_max), '\n']);

    fprintf(fid, ['alpha_UpperBoundaryCondition.type = BCT_CONSTANT_DERIVATIVE\n']);
    fprintf(fid, ['alpha_UpperBoundaryCondition.value = 0\n']);

    fprintf(fid, ['alpha_LowerBoundaryCondition.type = BCT_CONSTANT_DERIVATIVE\n']);
    fprintf(fid, ['alpha_LowerBoundaryCondition.value = 0\n']);
       
    fprintf(fid, ['L_UpperBoundaryCondition.type = BCT_CONSTANT_PSD\n']);
    fprintf(fid, ['L_UpperBoundaryCondition.value = 1e-21\n']);

    fprintf(fid, ['L_LowerBoundaryCondition.type = BCT_CONSTANT_VALUE\n']);
    fprintf(fid, ['L_LowerBoundaryCondition.value = 1e-21\n']);
    
    fprintf(fid, ['pc_LowerBoundaryCondition.type = BCT_CONSTANT_PSD\n']);
    fprintf(fid, ['pc_LowerBoundaryCondition.value = 1e-21\n']);

    fprintf(fid, ['pc_UpperBoundaryCondition.type = BCT_CONSTANT_DERIVATIVE\n']);
    fprintf(fid, ['pc_UpperBoundaryCondition.value = 1e-22\n']);
    
    fprintf(fid, ['psdRadialDiffusion.initial_PSD_Type = IPSDT_STEADY_STATE_L7\n']);
    fprintf(fid, ['psdRadialDiffusion.initial_PSD_Kp0 = 2\n']); %some base value to use throughout
    fprintf(fid, ['psdRadialDiffusion.initial_PSD_tauSteadyState = 2\n']);

    fclose(fid);

    system(['cat ', 'Parameters_DAonestep_Neumann.ini >> ',simulation_dir,'Parameters.ini']);

%% time and run VERB code
    cwd = pwd; %get current working directory
    code_time = tic;

    cd(simulation_dir); %move to simulation directory
    fprintf('Running VERB to compute SS initial condition...');
    system('./VERB_code_2_5_beta')
        
    fprintf('done.\n');
    cd(cwd); 

%% load PSD
    PSD = load_plt([simulation_dir,'Output/OutPSD.dat']);
    IPSD = PSD.arr(1,:,:,:);

%% Return to current working dir
    tstr = format_time(toc(code_time));
    fprintf('VERB RUN TIME: %s\n',tstr);
