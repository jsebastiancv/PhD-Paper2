%make sure it's marked
simulation_dir = [simulation_dir,'/'];
simulation_dir = regexprep(simulation_dir,'//','/');

%% check existance of simulation_dir directory
Vdir_check = exist(simulation_dir,'dir');
if Vdir_check ~= 7 %if does not exist then create it
    clear Vdir_check
    mkdir(simulation_dir);
    fprintf('%s',['copying files to ',simulation_dir,'...']);
    
    mkdir([simulation_dir,'Dxx_ini']);
    mkdir([simulation_dir,'Debug']);
    mkdir([simulation_dir,'Input']);
    mkdir([simulation_dir,'Output']);
    copyfile(['VERB/Dxx_ini/*'],[simulation_dir,'Dxx_ini/']);
    if copy_dxx_files %copy dxx ini and coeff files if need be
        Vdir_check2 = exist(source_simulation_dir,'dir');
        if Vdir_check2 ~= 7 %if does not exist then throw error
                error([source_simulation_dir,' does not exist! Can not copy dxx files\n']);
        else
            copyfile([source_simulation_dir,'Dxx_ini/*'],[simulation_dir,'Dxx_ini/']);
            try
                copyfile([source_simulation_dir,'DiffCoeff/*'],[simulation_dir,'DiffCoeff/']);          
            catch end
        end
    end
    fprintf('done\n')
else
    % let's check that we are not using the sim folder at the moment, and prompt the user to specify a new directory if necessary
    logchk = dir([simulation_dir,'logfile.log']);
    if ~isempty(logchk)
        dt = now - logchk.datenum;
        if dt < 10/86400
            errtxt = sprintf(['\n ******** SIM DIR ERROR ************\n',...
                ' logfile.log exists and was modified less than 10',...
                ' seconds ago.\n Are you trying to compute two simulations from the',...
                ' same directory?\n SIMULATION DIR is set to: %s.\n',...
                ' logfile.log modified: %s\n',...
                ' ***********************************\n'],simulation_dir,...
                now_time(logchk.datenum));
            error(errtxt)
        elseif dt < 5/1440 % if log file was modified less than 10 minutes ago then 
                       % prompt user to verify 
            proceed = false;
            while ~proceed
                wtxt = sprintf([' Simulation directory exists and logfile.log was',...
                    ' modified less than 5 minutes ago\n']);
                warning(wtxt);
                dtxt = sprintf(['Current simulation dir is: %s.\nDo you want to',...
                            ' continue using this directory? Y/N: '],simulation_dir);
                decision = lower(input(dtxt,'s'));
                if strcmp(decision,'n')
                    simulation_dir = input('Enter new simulation_dir name: ','s');
                    simulation_dir = [simulation_dir,'/'];
                    simulation_dir = regexprep(simulation_dir,'//','/');

                    logchk = dir([simulation_dir,'logfile.log']);
                    if ~isempty(logchk)
                        if dt < 10/86400
                            errtxt = sprintf(['\n ******** SIM DIR ERROR ************\n',...
                                ' logfile.log exists and was modified less than 10',...
                                ' seconds ago.\n Are you trying to compute two simulations from the',...
                                ' same directory?\n SIMULATION DIR is set to: %s.\n',...
                                ' logfile.log modified: %s\n',...
                                ' ***********************************\n'],simulation_dir,...
                                now_time(logchk.datenum));
                            error(errtxt)
                        elseif dt > 5/1440 % if log file was modified less than 10 minutes ago then 
                            proceed = true;
                        end
                    else
                        proceed = true;
                        Vdir_check = exist(simulation_dir,'dir');
                        if Vdir_check ~= 7 %if does not exist then create it
                            clear Vdir_check
                            mkdir(simulation_dir);
                            fprintf('%s',['copying files to ',simulation_dir,'...']);
                            mkdir([simulation_dir,'DiffCoeff']);
                            mkdir([simulation_dir,'Dxx_ini']);
                            mkdir([simulation_dir,'Debug']);
                            mkdir([simulation_dir,'Input']);
                            mkdir([simulation_dir,'Output']);
                            copyfile(['VERB/Dxx_ini/*'],[simulation_dir,'Dxx_ini/']);
                            if copy_dxx_files %copy dxx ini and coeff files if need be
                                Vdir_check2 = exist(source_simulation_dir,'dir');
                                if Vdir_check2 ~= 7 %if does not exist then throw error
                                        error([source_simulation_dir,' does not exist! Can not copy dxx files\n']);
                                else
                                    copyfile([source_simulation_dir,'Dxx_ini/*'],[simulation_dir,'Dxx_ini/']);
                                    copyfile([source_simulation_dir,'DiffCoeff/*'],[simulation_dir,'DiffCoeff/']);
                                end
                            end
                            fprintf('done\n')
                        end
                    end
                else
                    proceed = true;
                end
            end
        end
    end
end

copyfile(['VERB/VERB_code_2*'],[simulation_dir]);
