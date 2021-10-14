function [sw_data,sw_time] = get_omnimin_period(matlabd,varargin);

sw_data=[];
sw_time=[];

use_cache = false;
inx = 1:57;
if nargin > 1
    for it=1:1:nargin-1
        if strcmp(varargin(it),'use_cache')
            use_cache=true;
        elseif strcmp(varargin(it),'indices');
           inx = varargin{it+1};
           it = it+1;
        end
    end
end

[yr,~,~,~,~,~] = datevec(matlabd);

uniq_yr = unique(yr);

for iyear=1:length(uniq_yr);
    year = uniq_yr(iyear);
    if (use_cache)
        [sw_d,sw_t] = get_omnimin(year,'use_cache');
    else
        [sw_d,sw_t] = get_omnimin(year);
    end
    zz = find(sw_t >= matlabd(1)-1 & sw_t <= matlabd(end)+1);
    if ~isempty(zz)                
        sw_data = cat(1,sw_data,sw_d(zz,inx));
        sw_time = cat(1,sw_time,sw_t(zz));
    end
end 
sw_d = [];
for iinx=1:length(inx);
    if inx(iinx) == 47
        error('Kp is already interpolate, reload and interpolate to your specified time')
    else
        sw_d(:,iinx) = interp1(sw_time,sw_data(:,iinx),matlabd,'linear',NaN);
    end
end
sw_data = sw_d;
sw_time = matlabd;
