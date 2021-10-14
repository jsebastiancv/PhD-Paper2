%function used to load and process omni solar wind files
% Written by Adam Kellerman

function [sw_data, matlab_day] = get_omnimin(year,varargin)

global DAdata_dir SW_data_dir pll Satellite_data_dir

use_cache=false;
if nargin > 1
    for it=1:nargin-1
        if strcmp(varargin(it),'use_cache')
            use_cache=true;
            cache_dir=[getenv('HOME'),'/projects/cache/'];
            system(['mkdir -p ',cache_dir]);
        end
    end
end

ofdir = [Satellite_data_dir,'OMNI/processed_mat/'];
ofile = [num2str(year), '_onemin.mat']; 
matfile = dir([ofdir, ofile]); 
process = 0;
if ~isempty(matfile)
    %check date of omni.dat and omni.mat files, remake if dat is newer
    date_matfile = datenum(matfile.date);
    datfile = dir([Satellite_data_dir,'OMNI/omni_min', num2str(year), '.asc']);    

    date_datfile = datenum(datfile.date);
       
    if date_datfile > date_matfile
        fprintf('Dat file is newer: %s\n than mat file: %s\n',datfile.date,matfile.date);        
        process=1;
    else        
        fprintf('Dat file is older: %s\n than mat file: %s\n',datfile.date,matfile.date);
        fprintf('%s\n', ['Loading ',num2str(year), '_onemin.mat...']);
        if use_cache
            system(['rsync -avPzh ',ofdir,ofile,' ',cache_dir]);
            load([cache_dir,ofile]);
        else
            load([Satellite_data_dir,'OMNI/processed_mat/', num2str(year), '_onemin.mat']); 
        end        
        return
    end
else
    process=1;
    fprintf('%s\n', ['no file']);
end
if process 
    %%
    try
        parpool(4);
    catch
        fprintf('no parpool available or other error')
    end
    fprintf('%s', ['Searching for ','omni_min', num2str(year), '.asc','...']);
    txt_file = dir([Satellite_data_dir,'OMNI/omni_min', num2str(year), '.asc']);   
    if isempty(txt_file)
        fprintf('no asc file\n');
        sw_data = [];
        matlab_day = [];
    else
        fprintf('found\n');
        
        fprintf('%s', ['Loading ','omni_min', num2str(year), '.asc','...']);
        sw_data = load([Satellite_data_dir,'OMNI/omni_min', num2str(year), '.asc']);   
        sw_data_old = load([Satellite_data_dir,'OMNI/omni_min', num2str(year-1), '.asc']);    
        sw_data_old = sw_data_old(end-61:end, :);
        
        sw_data = cat(1, sw_data_old, sw_data);
        
        fprintf('%s\n', 'done.');
        %%
            
        fprintf('%s', 'Replacing fill values with NaN....');
        
        for i=1:46
            if (i == 4) || (i == 5) || (i == 6)
                z=find(sw_data(:,i) >= 99);
                sw_data(z,i) = NaN;
            end
            if (i == 7) || (i == 8) || (i == 9)
                z=find(sw_data(:,i) >= 999);
                sw_data(z,i) = NaN;
            end
            if (i == 10) || (i == 11) || (i == 13)
                z=find(sw_data(:,i) >= 999999);
                sw_data(z,i) = NaN;
            end
            if (i == 12) || (i == 28)
                z=find(sw_data(:,i) >= 99.9900);
                sw_data(z,i) = NaN;
            end
            if (i == 14) || (i == 15) || (i == 16) || (i == 17) ...
                    || (i == 18) || (i == 19) || (i == 20) || (i == 21) ...
                    || (i == 32) || (i == 33) || (i == 34) || (i == 35) ...
                    || (i == 36) || (i == 37)
                z=find(sw_data(:,i) >= 9999.99);
                sw_data(z,i) = NaN;
            end
            if (i == 22) || (i == 23) || (i == 24) || (i == 25) ...
                    z=find(sw_data(:,i) >= 99999.9);
                sw_data(z,i) = NaN;
            end
            if (i == 26) || (i == 29) || (i == 30) || (i == 31) || (i == 45)
                z=find(sw_data(:,i) >= 999.9900);
                sw_data(z,i) = NaN;
            end
            if (i == 27)
                z=find(sw_data(:,i) >= 9999999);
                sw_data(z,i) = NaN;
            end
            if (i == 38) || (i == 39) || (i == 40)
                z=find(sw_data(:,i) >= 999999);
                sw_data(z,i) = NaN;
            end
            if (i == 46)
                z=find(sw_data(:,i) >= 99.9000);
                sw_data(z,i) = NaN;
            end
            if (i == 42)
                z=find(sw_data(:,i) >= 99999);
                sw_data(z,i) = NaN;
            end
        end
        for j=1:46
            vtmp = sw_data(:,j);
            len = length(vtmp);
            xvals = 1:1:len;
            inan=find(isnan(vtmp));
            notnan=find(~isnan(vtmp));
            if ~isempty(inan) && length(notnan) > 1
                vtmp2 = interp1(xvals(notnan), vtmp(notnan), xvals, 'linear',NaN);
                sw_data(:,j) = vtmp2;
            end
            vtmp = sw_data(:,j);
            inan=find(isnan(vtmp));
            notnan=find(~isnan(vtmp));
        end
        fprintf('%s\n', 'done.');
        %%
        %convert sym-h to dst (Wanliss and Showalter, 2006)
        fprintf('%s', 'Calculating Dst from SYM-H (W&S 06)....');
        zz= find(sw_data(:,42) >= -300);
        A1 = 0.959;
        B1 = -2.5;
        A2 = 0.727;
        B2 = -62.2;
        
        % greater than -300 nT
        nel=length(sw_data(:,1));
        Dst = nan(1, nel);
        Dst(zz) = A1*sw_data(zz,42)+B1;
        % less than -300 nT
        zz= find(sw_data(:,42) < -300);
        Dst(zz) = A2*sw_data(zz,42)+B2;
        fprintf('%s\n', 'done.');

        %%
        matlab_day=datenum(sw_data(:,1),1,sw_data(:,2),sw_data(:,3),sw_data(:,4),0);
        %Kp = sw_data(:,47);
        Load_kp
        
        Kp = interp1(kp_time, kp_index, matlab_day, 'nearest'); %in actual format
        %%
        
        fprintf('%s\n', 'Calculating G and W....');    
        make_G_and_W

        sw_data = cat(2, sw_data, Kp);  %47
        sw_data = cat(2, sw_data, Dst');
        sw_data = cat(2, sw_data, G1);
        sw_data = cat(2, sw_data, G2);
        sw_data = cat(2, sw_data, G3);
        sw_data = cat(2, sw_data, W1);
        sw_data = cat(2, sw_data, W2);
        sw_data = cat(2, sw_data, W3);
        sw_data = cat(2, sw_data, W4);
        sw_data = cat(2, sw_data, W5);
        sw_data = cat(2, sw_data, W6); %57
        
        fprintf('%s\n', 'done');
        %%
        fprintf('%s', 'Resizing array to year only..'); 
        zz=find(sw_data(:,1) == year);
        sw_data = sw_data(zz,:);
        matlab_day = matlab_day(zz);
        fprintf('%s\n', 'Done...');
        %%
        fprintf('%s', 'Saving....'); 
    %    save([DAdata_dir, 'Processed_data/omni_data/', num2str(year), '_onemin.mat'], 'sw_data', 'matlab_day')    
    %   save(['/home/akellerman/Documents/omni_data/', num2str(year), '_onemin.mat'], 'sw_data', 'matlab_day')    
        save([Satellite_data_dir,'OMNI/processed_mat/', num2str(year), '_onemin.mat'], 'sw_data', 'matlab_day');
        fprintf('%s\n', 'done.');
    end
end
