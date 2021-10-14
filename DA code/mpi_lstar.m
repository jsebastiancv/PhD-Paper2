function [Flux, Flux_at_R] = run_mpi(SimTime,SimPSD,SimL,SimInvMu,SimInvK,target_MLT, target_epc,...
    target_alpha, mpi_dir, lt_date, kext,numcpu,recalc_mpi);

recalc = false;
ntall = length(SimTime);
target_pc = pfunc(target_epc);
datetag = ['test']; %traditionally used to specify the date, though you can put whatever you want here.
sate = 'Reanalysis';
rvals = 1:0.25:10;
n_pos = length(rvals);

if kext == 4
    mfm = 'T89';
elseif kext == 11
    mfm = 'T04s';
elseif kext == 13
    mfm = 'TS07D';
end

file = dir([mpi_dir,'/output/',sate,'_LstarandK_',mfm,'_',datetag,'.txt']);
file2 = dir([mpi_dir,'/output/',sate,'_Bdata_',mfm,'_',datetag,'.txt']);
if ~isempty(file2) && ~recalc_mpi
    if file2.datenum < lt_date
        recalc = true;
        pause
    else
        lstar_data = load([mpi_dir,'/output/',file.name]);
        alpha_eq = reshape(lstar_data(:,1),n_pos,[])';
        B_data = load([mpi_dir,'/output/',file2.name]);
        B_time = reshape(B_data(:,1),n_pos,[])';
        B_time = B_time(:,1) ./ 86400 + datenum(1970,1,1);
        B_MLT = reshape(B_data(:,3),n_pos,[])';
        ttol = 1/86400; % tolerance for time, it's usually 0, but just in case
        atol = 0.001; % tolerance for alpha_eq, same as above
        mtol = 0.1; % mlt tolerance, this one has an error ~0.01 in most cases

        % dealing with mean of MLT.
        B_MLT = reshape(B_MLT,[],1);
        x=cos(B_MLT./24 .* 2.*pi);
        y=sin(B_MLT./24 .* 2.*pi);
        ang = atan(mean(y)/mean(x));
        if ang < 0
            ang = ang+2.*pi;
        end
        B_MLT = ang.*24 ./ (2.*pi);

        %abs(B_time(1) - SimTime(1))
        %abs(B_time(end) - SimTime(end))
        %length(B_time)
        %abs(alpha_eq(1) - target_alpha)
        %target_MLT
        %abs(B_MLT - target_MLT)
        if abs(B_time(1) - SimTime(1)) > ttol || ...
            abs(B_time(end) - SimTime(end)) > ttol || ...
            length(B_time) ~= ntall || ...
            abs(alpha_eq(1) - target_alpha) > atol || ...
            abs(B_MLT - target_MLT) > mtol
            recalc = true;
        end
    end
end
if recalc_mpi
    recalc = true;
end
if isempty(file) || recalc
    sysaxes = 4; 
    if strcmp(mfm,'TS07d')
        fprintf('WARNING!! the TS07d model may take a very long time to compute even one time point\n') 
    end
    ang = target_MLT./24 .*2 .*pi +pi; % starting from 12 as x and y are positive 
        % from there.   Note that this is a rough approximation, we could actually
        % use the MLT from the model in the future
    if ang > 2.*pi
        ang = ang - (2.*pi);
    end
    xvals = rvals .* cos(ang);
    yvals = rvals .* sin(ang);
    zvals = xvals.*0; 


    xSM = cat(1,xvals,yvals,zvals)';

    maginput = make_maginput_all(SimTime); % call function to make inputs
    maginput = maginput';

    lsw_dat=nan(ntall,29);

    tvec = datevec(SimTime);
    doy  = date2doy(SimTime);

    lsw_dat(:,1) = tvec(:,1);
    lsw_dat(:,2) = doy;
    lsw_dat(:,3:7) = tvec(:,2:6);
    lsw_dat(:,8) = (SimTime-datenum(1970,1,1)).*86400;
    %lsw_dat(:,9:11) = [0 0 0] % position comes from a file now
    lsw_dat(:,12:29) = NaN;
    lsw_dat(:,19) = maginput(:,1);
    lsw_dat(:,18) = maginput(:,2);
    lsw_dat(:,16) = maginput(:,3);
    lsw_dat(:,15) = maginput(:,4);
    lsw_dat(:,17) = maginput(:,5);
    lsw_dat(:,13) = maginput(:,6);
    lsw_dat(:,14) = maginput(:,7);
    lsw_dat(:,21:29) = maginput(:,8:16);

    for isw=12:29
        lsw_dat(isnan(lsw_dat(:,isw)),isw) = -1.0d31;
    end

    lsw_dat(isnan(lsw_dat))=-1.0d31;

    svdir=[mpi_dir,'/input/'];
    save([svdir,'target_alpha_',datetag,'.txt'],'target_alpha','-ASCII','-DOUBLE');
    save([svdir,'SM_positions_',datetag,'.txt'],'xSM','-ASCII','-DOUBLE');

    file_out = [sate,'_SM_LocationandSolarWindData_',datetag,'.txt'];
    save([mpi_dir,'/input/', file_out],'lsw_dat','-ASCII','-double');
    [~,~]=system(['gzip -f --fast ',mpi_dir,'/input/',file_out]);

    str_kext = num2str(kext);
    str_sysaxes = num2str(sysaxes);
    str_numcpu = num2str(numcpu);
    mfck = dir([mpi_dir,'/mpi_lstar']);
    if isempty(mfck)
        cwd = pwd;
        cd(mpi_dir)
        % try to compile the executible
        system('mpif90 mpi_lstar_BLK.f90 LstarCalc.f90 -Llib -loneradesp_gnu64_ntmax1 -o mpi_lstar -O3');
        mfck = dir([mpi_dir,'/mpi_lstar']);
        if isempty(mfck)
            cd(cwd)
            error('MPI compile failed, check/recompile library file IRBEM, name it appropriately, and put in lib/')
        end
        cd(cwd)
    end    
    system([mpi_dir,'/run_mpi.sh ',sate,' ',datetag,' ',str_kext,' ',...
        str_sysaxes,' ',str_numcpu]);

    files = dir([mpi_dir,'/output/',sate,'_LstarandK_',mfm,'_',datetag,'.txt.gz']);
    files2 = dir([mpi_dir,'/output/',sate,'_Bdata_',mfm,'_',datetag,'.txt.gz']);
    system(['gunzip -vf ',mpi_dir,'/output/',files.name]);
    system(['gunzip -vf ',mpi_dir,'/output/',files2.name]);

    file = dir([mpi_dir,'/output/',sate,'_LstarandK_',mfm,'_',datetag,'.txt']);
    file2 = dir([mpi_dir,'/output/',sate,'_Bdata_',mfm,'_',datetag,'.txt']);
end

lstar_data = load([mpi_dir,'/output/',file.name]);
lstar_data(lstar_data < -1e29) = NaN;
Lstar = reshape(lstar_data(:,2),n_pos,[])';
InvK = reshape(lstar_data(:,3),n_pos,[])';
alpha_eq = reshape(lstar_data(:,1),n_pos,[])';% should be same as 1

B_data = load([mpi_dir,'/output/',file2.name]);
B_data(B_data < -1e29) = NaN;

B_equator = reshape(B_data(:,4),n_pos,[])'; %field on the equator (found in FORTRAN) 
B_equator = B_equator * 1e-5;
MLT = reshape(B_data(:,3),n_pos,[])';

Lstar(Lstar <= 0) = NaN;
InvK(InvK < 0) = NaN;
B_equator(B_equator <= 0) = NaN;
alpha_eq(alpha_eq < 0) = NaN;

mc2 = 0.511;
InvMu = nan(ntall,size(Lstar,2));
Flux = nan(ntall,size(Lstar,2));
nl = size(Lstar,2);
for il=1:nl
    InvMu(:,il) = (target_pc .* ...
        sin(target_alpha .* pi/180)) ...
        .^2 ./ (B_equator(:,il) .* 2 * mc2);
end

SimL_vec = SimL(:,1,1);

vecSimInvMu = log10(reshape(SimInvMu,[],1));
vecSimInvK = reshape(SimInvK,[],1).^(1/3);
vecSimL = reshape(SimL,[],1);

szl = size(Lstar,2);
Flux_at_R = nan(ntall,szl);
Flux = nan(ntall,length(SimL_vec));
try
    parpool(numcpu);
catch end

%% interpolation part
parfor it = 1: ntall
    loop_counter(it,size(SimPSD,1),'tput','pretext','interpolation: ');
    fix_InvMu = reshape(InvMu(it,:),[],1);
    fix_InvK = reshape(InvK(it,:),[],1);
    fix_L = reshape(Lstar(it,:),[],1);
    Fluxt = nan(1,szl);
    zz=find(~isnan(fix_L) & ~isnan(fix_InvMu) & ~isnan(fix_InvK) & fix_L >= 1);
    if length(zz) > 1
        unql = unique(fix_L(zz));
        if length(unql) > 1
            R_PSD = log10(reshape(SimPSD(it,:,:,:),[],1));
            fix_PSD = 10.^(griddata(vecSimInvMu,vecSimInvK,vecSimL,R_PSD,...
                log10(fix_InvMu(zz)), fix_InvK(zz).^(1/3),...
                fix_L(zz),'natural'));
            Fluxt(zz) = fix_PSD .* target_pc.^2;
            Flux(it,:) = interp1(fix_L(zz),Fluxt(zz),SimL_vec,'linear');
        end
    end
    Flux_at_R(it,:) = Fluxt;
end

Flux_at_R(Flux_at_R <= 0) = 1e-22;
Flux(Flux <= 0) = 1e-22;

