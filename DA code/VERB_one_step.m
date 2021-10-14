%% save initial PSD for VERB computation
    init_PSD_file = [simulation_dir,'Input/da3d_init_psd.dat'];
    save_plt_lpa(init_PSD_file, lastPSD);

%% parameters structure:
    steps = 1;
    fid = fopen([simulation_dir,'Parameters.ini'], 'w');
    fprintf(fid, ['timeStep = ', num2str(timestep), '\n']);
    
    fprintf(fid, ['nDays = ', num2str(d_time), '\n']);
    fprintf(fid, ['general_Output_parameters.iterStep = ',num2str(steps), '\n']);
        
    fprintf(fid, 'useKp = constant\n');
    fprintf(fid, ['constKp = ', num2str(Kp_sim), '\n']);
    
%     fprintf(fid, 'useBf = constant\n');
%     fprintf(fid, ['constBf = ', num2str(Bf_sim_rec), '\n']);
    
    fprintf(fid, 'useBf =  No\n');
        
    fprintf(fid, 'useLpp = constant\n');
    fprintf(fid, ['constLpp = ', num2str(Lpp_sim), '\n']);
    
    fprintf(fid, 'useRadialDiffusion = Yes\n');    
    fprintf(fid, 'usePcAlphaMixedTerms = Yes \n');
    fprintf(fid, 'useEnergyDiffusion = Yes \n');
    fprintf(fid, 'useAlphaDiffusion = Yes \n');
    
    fprintf(fid, ['localDiffusionsGrid_L.size = ', num2str(NL), '\n']);
    fprintf(fid, ['localDiffusionsGrid_L.min = ', num2str(L_min), '\n']);
    fprintf(fid, ['localDiffusionsGrid_L.max = ', num2str(L_max), '\n']);

    fprintf(fid, ['localDiffusionsGrid_epc.size = ', num2str(NE), '\n']);
    fprintf(fid, ['localDiffusionsGrid_epc.min = ', num2str(E_min), '\n']);
    fprintf(fid, ['localDiffusionsGrid_epc.max = ', num2str(E_max), '\n']);

    fprintf(fid, ['localDiffusionsGrid_alpha.size = ', num2str(NA), '\n']);
    fprintf(fid, ['localDiffusionsGrid_alpha.min = ', num2str(a_min), '\n']);
    fprintf(fid, ['localDiffusionsGrid_alpha.max = ', num2str(a_max), '\n']);
    
%     if mp_set % mp losses, lossy outer boundary
%         fprintf(fid, ['L_UpperBoundaryCondition.type =  BCT_CONSTANT_VALUE\n']);
%         fprintf(fid, ['L_UpperBoundaryCondition.value = 1e-21\n']);        
%     else
%         fprintf(fid, ['L_UpperBoundaryCondition.type =  BCT_CONSTANT_DERIVATIVE\n']);
%         fprintf(fid, ['L_UpperBoundaryCondition.value = 0\n']);
%     end
    
    fprintf(fid, ['alpha_UpperBoundaryCondition.type = BCT_CONSTANT_DERIVATIVE\n']);
    fprintf(fid, ['alpha_UpperBoundaryCondition.value = 0\n']);

    fprintf(fid, ['alpha_LowerBoundaryCondition.type = BCT_CONSTANT_DERIVATIVE\n']);
    fprintf(fid, ['alpha_LowerBoundaryCondition.value = 0\n']);
%        
%     fprintf(fid, ['L_UpperBoundaryCondition.type = BCT_CONSTANT_PSD\n']);
%     fprintf(fid, ['L_UpperBoundaryCondition.value = 1e-21\n']);

    fprintf(fid, ['L_LowerBoundaryCondition.type = BCT_CONSTANT_VALUE\n']);
    fprintf(fid, ['L_LowerBoundaryCondition.value = 1e-21\n']);

    fprintf(fid, ['L_UpperBoundaryCondition.type = BCT_FILE_GRID\n']);
    fprintf(fid, ['L_UpperBoundaryCondition.filename = ./Input/L_upper_BC.dat\n']);
    
    fprintf(fid, ['pc_LowerBoundaryCondition.type = BCT_CONSTANT_PSD\n']);
    fprintf(fid, ['pc_LowerBoundaryCondition.value = 1e-21\n']);

    fprintf(fid, ['pc_UpperBoundaryCondition.type = BCT_CONSTANT_DERIVATIVE\n']);
    fprintf(fid, ['pc_UpperBoundaryCondition.value = 1e-22\n']);

    fprintf(fid, ['psdRadialDiffusion.initial_PSD_Type = IPSDT_FILE\n']);
    fprintf(fid, ['psdRadialDiffusion.initial_PSD_fileName = ./Input/da3d_init_psd.dat\n']);
   
    fprintf(fid, ['sourcesAndLosses.useLmp = true\n']);
     fprintf(fid, ['sourcesAndLosses.fileLmp = ./mp_sim.txt\n']);
%     
    if Pdyn_sim > 3
        fprintf(fid, ['CoeffFileName = ./Dxx_ini_2016/Daa_EMIC_plume.ini\n']);
%         fprintf(fid, ['CoeffFileName = ./Dxx_ini_2016/Daa_EMIC_plume_pp.ini\n']);
    end   
%     
    system(['cat ', 'Parameters_DAonestep_Neumann.ini >> ',simulation_dir,'Parameters.ini']);

%% time and run VERB code

    cwd = pwd; %get current working directory
    code_time = tic;
    
    cd(simulation_dir);

     fileID = fopen('mp_sim.txt','w');
     fprintf(fileID,'%e   ',Mp_it);
     fclose(fileID);
%     
     fileID = fopen('mp_sim.txt','a');
     fprintf(fileID,'\n');
     fprintf(fileID,'%e   ',extra);
     fclose(fileID);

    fprintf('Running one-step VERB ...');
%     system('./VERB_code_2_5_beta')
    system('./VERB_code_2_4_1_DA_2020_newDxx')
    fprintf('done.\n');
    
    cd(cwd); 

%% Return to current working dir
    tstr = format_time(toc(code_time));
    fprintf('VERB RUN TIME: %s\n',tstr);
