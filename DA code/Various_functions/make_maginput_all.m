% written by A. Kellerman, 2014

function [maginput] = make_maginput_all(matlabd,varargin)

% uses both realtime and historic data, depending on availability
global Satellite_data_dir

use_cache=false;
if nargin > 1
    for it=1:nargin-1
        if strcmp(varargin(it),'use_cache')
            use_cache=true;
        end
    end
end
nT = length(matlabd);

st = matlabd(1);
et = matlabd(end);

[syr,~,~,~,~,~] = datevec(st);
[eyr,~,~,~,~,~] = datevec(et);

sw_d=[];
sw_t=[];
if syr ~= eyr
    for iyr = syr:eyr
        if (use_cache)
            [sw_data, sw_matlabd]=get_omnimin(iyr,'use_cache');
        else
            [sw_data, sw_matlabd]=get_omnimin(iyr);
        end
        sw_d = cat(1,sw_d,sw_data);
        sw_t = cat(1,sw_t,sw_matlabd);
    end
    sw_data = sw_d;
    sw_matlabd = sw_t;
    clear sw_d sw_t
else
    if (use_cache)
        [sw_data, sw_matlabd]=get_omnimin(syr,'use_cache');
    else
        [sw_data, sw_matlabd]=get_omnimin(syr);
    end
end


%% arrays are different for omni and realtime so we combine here
if ~isempty(sw_data) % empty if asc file does not exist
    zz = find(~isnan(sw_data(:,22)) & ~isnan(sw_data(:,26)) & ~isnan(sw_data(:,19)) ...
    & ~isnan(sw_data(:,28)));

    sw_data = sw_data(zz,:);
    sw_matlabd = sw_matlabd(zz,:);

    
    try
        load([Satellite_data_dir,'/ACE/sw_data_prop_simple.mat']);
        sw_data_rt = shifted_sw_data;
        clear shifted_sw_data % 3-5 B GSM, 7 den, 8 vel, 9 Pdyn

        zz=find(sw_data_rt(:,1) > sw_matlabd(end));
        sw_data_rt = sw_data_rt(zz,:);
        nanpad = nan(size(sw_data_rt,1),size(sw_data,2));
        sinx = size(sw_data,1) + 1;
        sw_data = cat(1,sw_data,nanpad); % pad with nans then replace below

        sw_matlabd = cat(1,sw_matlabd,sw_data_rt(:,1));
        sw_data(sinx:end,15) = sw_data_rt(:,3);
        sw_data(sinx:end,18) = sw_data_rt(:,4);
        sw_data(sinx:end,19) = sw_data_rt(:,5);
        sw_data(sinx:end,26) = sw_data_rt(:,7);
        sw_data(sinx:end,22) = sw_data_rt(:,8);
        sw_data(sinx:end,28) = sw_data_rt(:,9);
    catch
        fprintf('no realtime propagated plasma available, or path is wrong\n');
    end
else
    try
        load([Satellite_data_dir,'/ACE/sw_data_prop_simple.mat']);
        sw_data_rt = shifted_sw_data;
        clear shifted_sw_data % 3-5 B GSM, 7 den, 8 vel, 9 Pdyn
        sw_data = nan(size(sw_data_rt,1),57);
        sw_matlabd = sw_data_rt(:,1);
        sw_data(:,15) = sw_data_rt(:,3); % imf x
        sw_data(:,18) = sw_data_rt(:,4); % imf y
        sw_data(:,19) = sw_data_rt(:,5); % imf z
        sw_data(:,26) = sw_data_rt(:,7); % Den
        sw_data(:,22) = sw_data_rt(:,8); % Vel
        sw_data(:,28) = sw_data_rt(:,9); % Pdyn
    catch
        fprintf('Realtime propagated plasma data path wrong, skipping\n');
    end
end
try
    load([Satellite_data_dir,'/Dst/Dst.mat']);
    Dst_time = Dst(:,1);
    Dst_current = Dst(:,2);
    Dst_forecast = Dst(:,3);
    zz=find(isnan(sw_data(:,48)));
    sw_data(zz,48) = interp1(Dst_time,Dst_current,sw_matlabd(zz),'linear',NaN);
catch
    fprintf('Realtime Dst file path is wrong, skipping\n');
end

zz=find(sw_matlabd >= st - 1 & sw_matlabd <= et + 1);
sw_data = sw_data(zz,:);
sw_matlabd = sw_matlabd(zz);

make_G_and_W %both omni and rt are in minute resolution

zz= find(isnan(sw_data(:,49)));
if ~isempty(zz)
    sw_data(zz,49) = G1(zz);
end
zz= find(isnan(sw_data(:,50)));
if ~isempty(zz)
    sw_data(zz,50) = G2(zz);
end
zz= find(isnan(sw_data(:,51)));
if ~isempty(zz)
    sw_data(zz,51) = G3(zz);
end
for iw=1:6
    num=num2str(iw+51);
    zz= find(isnan(sw_data(:,iw+51)));
    if ~isempty(zz)
        eval(['sw_data(zz,',num,') = W',num2str(iw),'(zz);']);
    end
end
maginput=nan(25,nT);

%kp
Load_kp
kp_index_f = kp_index;
kp_time_f = kp_time;
if (true)

    %% let's try to load Wing Kp
    try
        kpdata = importdata([Satellite_data_dir,...
                'Kp/Wing_Predicted_Kp/wing-kp.txt'],' ',20);
        inx = 1 ; 
        Wing_Kp_time1 = nan(length(kpdata.data),1);
        Wing_Kp_index1 = nan(length(kpdata.data),1);
        Wing_Kp_time4 = nan(length(kpdata.data),1);
        Wing_Kp_index4 = nan(length(kpdata.data),1);
        for it = 1:length(kpdata.data)
            hhmm = sprintf('%0.4i\n',kpdata.data(it,8));
            Wing_Kp_time1(inx) = datenum(kpdata.data(it,5),kpdata.data(it,6),kpdata.data(it,7),...
                str2num(hhmm(1:2)),str2num(hhmm(3:4)),0);
            Wing_Kp_index1(inx) = round(kpdata.data(it,9) .* 10);

            hhmm = sprintf('%0.4i\n',kpdata.data(it,13));
            Wing_Kp_time4(inx) = datenum(kpdata.data(it,10),kpdata.data(it,11),kpdata.data(it,12),...
                str2num(hhmm(1:2)),str2num(hhmm(3:4)),0);
            Wing_Kp_index4(inx) = round(kpdata.data(it,14) .* 10);
            inx=inx+1;
        end    
        if Wing_Kp_time1(1) > kp_time_f(end)    
            tmpkp_time = kp_time_f(end):3/24:Wing_Kp_time1(1);
            tmpkp_time = tmpkp_time(2:end-1)';
            tmpkp = tmpkp_time+NaN;
            kp_index_f = cat(1,kp_index_f,tmpkp);
            kp_time_f = cat(1,kp_time_f,tmpkp_time);
        end

        zz=find(Wing_Kp_time1 > kp_time_f(end));
        kp_index_f = cat(1,kp_index_f,Wing_Kp_index1(zz));
        kp_time_f = cat(1,kp_time_f,Wing_Kp_time1(zz));

        zz=find(Wing_Kp_time4 > kp_time_f(end));
        kp_index_f = cat(1,kp_index_f,Wing_Kp_index4(zz));
        kp_time_f = cat(1,kp_time_f,Wing_Kp_time4(zz));
    catch
        fprintf('No Wing Kp included\n');
    end
end
try
    load([Satellite_data_dir,'Kp/Kp_3day/kp_3day.mat']);
    zz = find(isnan(kp_index_f));
    if ~isempty(zz)
        kp_index_f(zz) = interp1(kp_time,kp_index.*10,kp_time_f(zz),'nearest',NaN);
    end    
    if kp_time(1) > kp_time_f(end)    
        tmpkp_time = kp_time_f(end):3/24:kp_time(1);
        tmpkp_time = tmpkp_time(2:end-1)';
        tmpkp = tmpkp_time+NaN;
        kp_index = cat(1,kp_index_f,tmpkp,kp_index.*10);
        kp_time = cat(1,kp_time_f,tmpkp_time,kp_time);
    else
        zz=find(kp_time > kp_time_f(end));
        kp_index = cat(1,kp_index_f,kp_index(zz).*10);
        kp_time = cat(1,kp_time_f,kp_time(zz));
    end
catch
    fprintf('no 3-day forecast kp available\n');
end

Kp = interp1(kp_time,kp_index,matlabd,'nearest',NaN);
maginput(1,:) = Kp;

dt = median(diff(matlabd))/2.;

tmplsw = zeros(7,nT)+NaN;
parfor it=1:nT
    loop_counter(it,nT,'tput','mod',1000,'pretext','maginput, mag plasma: ');
    svar=0;
    zz=find(sw_matlabd >= matlabd(it)-dt & sw_matlabd < matlabd(it)+dt);
    if ~isempty(zz)
        tmpvars = zeros(7,1)+NaN;
        for ivar=1:7

            if ivar == 7
                svar=15; %Bx, GSM

            elseif ivar == 5
                svar=18; %By, GSM

            elseif ivar == 6
                svar=19; %Bz, GSM

            elseif ivar == 3
                svar=22; %Flow speed, km/s

            elseif ivar == 2
                svar=26; %Proton Density, n/cc

            elseif ivar == 4
                svar=28; %Flow Pressure, nPa

            elseif ivar == 1
                svar=48; %Sym-H converted to Dst / Forecast Dst
            end
    
            tmpvals = sw_data(zz,svar);
            zz2=find(~isnan(tmpvals) & tmpvals ~=1e-32);
            if ~isempty(zz2)
                %lsw_dat(it,ivar+11) = mean(tmpvals(zz2));
                tmpvars(ivar) = mean(tmpvals(zz2));
            end
        end
        tmplsw(:,it) = tmpvars;
    end
end
maginput(2:7,:) = tmplsw(1:6,:);
clear tmplsw

tmplsw = zeros(1,nT);
parfor it=1:nT
    loop_counter(it,nT,'tput','mod',1000,'pretext','maginput,AL: ');
    zz=find(sw_matlabd >= matlabd(it)-dt & sw_matlabd < matlabd(it)+dt);
    if ~isempty(zz)
        tmpvals = sw_data(zz,39); %AL
        zz2=find(~isnan(tmpvals) & tmpvals ~=1e-32);
        if ~isempty(zz2)
            tmplsw(it) = mean(tmpvals(zz2));
        end
    end
end
maginput(17,:) = tmplsw;
clear tmplsw

tmplsw = zeros(9,nT)+NaN;
parfor it=1:nT
    loop_counter(it,nT,'tput','mod',1000,'pretext','maginput, Dst..: ');
    svar=0
    zz=find(sw_matlabd >= matlabd(it)-dt & sw_matlabd < matlabd(it)+dt);
    if ~isempty(zz)
        tmpvars=zeros(9,1)+NaN;
        for ivar=1:9
            svar=48+ivar;
            tmpval = sw_data(zz,svar);
            zz2=find(~isnan(tmpval) & tmpval ~=1e-32);
            if ~isempty(zz2)
                %                   lsw_dat(it,ivar) = mean(tmpval(zz2));
                tmpvars(ivar) = mean(tmpval(zz2));
            end
        end
        tmplsw(:,it) = tmpvars;
    end
end
maginput(8:16,:) = tmplsw;
maginput=double(maginput);
clear tmplsw
