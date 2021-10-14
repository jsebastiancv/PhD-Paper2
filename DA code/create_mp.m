function status = create_mp(dir,sdate,edate)
%% Creates mp parameter file for event
%

import plasma.Dip.Lalpha2K

%path2data = '/data/MP/MP/Processed_Mat_Files/';
path2data = '/export/mag5/rbm/data/rbm-data/MP/MP/Processed_Mat_Files/';
%mfm = 'TS07Dmid15';
mfm = 'T04s';
version = 'ver4';

dt = 1 / 24; % VERB time step, days
pa.size = 91; % size of pitch angle grid
pa.min = 0.7 * pi / 180; % minimum pitch angle, rad
pa.max = 89.3 * pi / 180; % maximum pitch angle, rad
pa.Lval = 6.6; % value at which K-grid is calculated in the VERB code
pa.use_log = true;

use_eq_lcds = false; % use 90-deg LCDS for all pitchangles

file_fstr = 'MP_nopos_LCDS2_%sto%s_%s_%s_%s.mat';

%% Read and process LCDS data
% if ~use_eq_lcds
%     error('The option has not yet been implemented');
% end

syear = year(datetime(datestr(sdate)));
smonth = month(datetime(datestr(sdate)));

eyear = year(datetime(datestr(edate)));
emonth = month(datetime(datestr(edate)));
eday = eomday(eyear, emonth);

sdate2 = datenum(syear, smonth, 1);
edate2 = datenum(eyear, emonth, eday);

cdate = sdate2;
lcds_all.arr = [];
lcds_all.time = [];
kinv_all.arr = [];
kinv_all.time = [];
while cdate < edate2
    cdate_next = addtodate(cdate, 1, 'month');
    
    current_file_lstar = sprintf(file_fstr, datestr(cdate, 'yyyymmdd'), ...
        datestr(cdate_next-1, 'yyyymmdd'), 'lstar', mfm, version);
    disp(['Reading file: ', current_file_lstar]);
    load([path2data, current_file_lstar]);
    lcds_all.arr = cat(1, lcds_all.arr, Lstar);
    lcds_all.time = cat(1, lcds_all.time, time);
        
    current_file_Kinv = sprintf(file_fstr, datestr(cdate, 'yyyymmdd'), ...
        datestr(cdate_next-1, 'yyyymmdd'), 'invk', mfm, version);
    disp(['Reading file: ', current_file_Kinv]);
    load([path2data, current_file_Kinv]);
    kinv_all.arr = cat(1, kinv_all.arr, InvK);
    kinv_all.time = cat(1, kinv_all.time, time);

    cdate = cdate_next;
end

if any(lcds_all.time ~= kinv_all.time)
    error('Lstar and K times are not consistent!');
end

% idx = (lstar_all.time >= sdate) & (lstar_all.time <= edate);
% lstar_all.arrr = lstar_all.arr(idx,:);
% lstar_all.time = lstar_all.time(idx);

time_verb = [sdate:dt:edate]';
nt_verb = length(time_verb);
if use_eq_lcds 
    disp('Using equatorial LCDS for all pitch angles');
    disp('Binning LCDS...');
    lcds_s.arr = squeeze(lcds_all.arr(:, end));
    lcds_s.time = lcds_all.time;
    
    lcds_t.time = time_verb;
    lcds_t.arr = nan(nt_verb, 1);
    for it = 1:nt_verb
        current_time = lcds_t.time(it);
        idx = (lcds_s.time >= current_time - dt/2) & ...
            (lcds_s.time <= current_time + dt/2);
        lcds_t.arr(it) = nanmean(lcds_s.arr(idx));
    end
    
    lcds_out = repmat(lcds_t.arr, [1, pa.size]);
    lcds_out = cat(2, lcds_t.time - sdate, lcds_out);
    disp('Done');
else
    disp('Interpolating LCDS into the VERB K-grid.');
    if pa.use_log
        pa.grid = 10 .^ linspace(log10(pa.min), log10(pa.max), pa.size);
    else
        pa.grid = linspace(pa.min, pa.max, pa.size);
    end
    pa.K_grid = Lalpha2K(pa.Lval, pa.grid);
        
    nt_obs = length(lcds_all.time);
    lcds_interp.time = lcds_all.time;
    lcds_interp.arr = nan(nt_obs, pa.size);
    for it = 1:nt_obs
        notnan = ~isnan(lcds_all.arr(it,:));
        if sum(notnan(:)) < 2
            continue;
        end
        lcds_interp.arr(it,:) = interp1(kinv_all.arr(it,notnan), ...
            lcds_all.arr(it,notnan), pa.K_grid, 'linear', NaN);
        nanidx = isnan(lcds_interp.arr(it,:));
        lcds_interp.arr(it,nanidx) = interp1(pa.K_grid(~nanidx), ...
            lcds_interp.arr(it,~nanidx), pa.K_grid(nanidx), 'nearest', 'extrap');
    end
    if any(isnan(lcds_interp.arr(:)))
        warning('There are NaNs in interpolated LCDS. Please extrapolate over them.');
    end
    
    lcds_t.time = time_verb;
    lcds_t.arr = nan(nt_verb, pa.size);
    for it = 1:nt_verb
        current_time = lcds_t.time(it);
        idx = (lcds_interp.time >= current_time - dt/2) & ...
            (lcds_interp.time <= current_time + dt/2);
        lcds_t.arr(it,:) = squeeze(nanmean(lcds_interp.arr(idx,:), 1));
    end
    
    lcds_out = cat(2, lcds_t.time - sdate, lcds_t.arr);
    disp('Done.');
end

if any(isnan(lcds_out(:)))
    error('There are NaNs in LCDS array. Please modify the code.');
end

save([dir 'Input/' 'mp.txt'], 'lcds_out', '-ascii');

status = 1;

end
