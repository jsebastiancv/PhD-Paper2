function [varargout]=load_plt(filename, varargin)

nozone = false;
dosqueeze = true;
dopermute = true;
first_zone = 1;
n_zones = 4096;
skip_zones = 0;

comments={};

if (nargin > 1)
    for it = 1:1:nargin-1
        if strcmp(varargin(it), 'nopermute')
            dopermute = false;
        elseif strcmp(varargin(it), 'permute')
            dopermute = true;
%            permute_mask = varargin(it+1);
        elseif strcmp(varargin(it), 'squeeze')
            dosqueeze = true;
        elseif strcmp(varargin(it), 'nosqueeze')
            dosqueeze = false;
        elseif strcmp(varargin(it), 'nozone')
            nozone = true;
        elseif strcmp(varargin(it), 'first_zone')
            first_zone = varargin{it+1};
            it = it + 1;
        elseif strcmp(varargin(it), 'n_zones')
            n_zones = varargin{it+1};
            it = it + 1;
        elseif strcmp(varargin(it), 'skip_zones')
            skip_zones = varargin{it+1};
            it = it + 1;
        end
    end
end
%    if (~(strcmp(varargin(1), 'dosqueeze')))
%        dosqueeze = false;
%    end;
%end;


fid = fopen(filename);

f_names = fgetl(fid);
while (length(f_names) == 0 || f_names(1) == '#')
    comments = [comments, f_names];
    f_names = fgetl(fid);
end;

f_names = regexpi(f_names, '"(.*?)"', 'match');
n_of_vars = size(f_names,2);
for it = 1:n_of_vars
%%%    vars(it).name = strrep(f_names(it){1}, '"', ''); %%%
    vars(it).name = strrep(f_names{it}, '"', ''); %%%
    vars(it).comments = comments;
end

zone_number = 0;
sz = 0;
zn = 0;
time = 0;
tmp = fgetl(fid);
while (~nozone && ((length(tmp) < 4) || (strcmpi(tmp([1:4]), 'ZONE') ~= 1)))
    tmp = fgetl(fid);
end;

% assume all zones have the same size. Would not work otherwise anyway.
t = regexpi(tmp, 'T="(.*)"', 'tokens');      
if (length(t) > 0) time = str2num(t{1}{1}); end;

t = regexpi(tmp, 'I=\D*([0-9]+)', 'tokens');
if (length(t) > 0) size1 = str2num(t{1}{1}); else size1 = 0; end;
t = regexpi(tmp, 'J=\D*([0-9]+)', 'tokens');
if (length(t) > 0) size2 = str2num(t{1}{1}); else size2 = 0; end;       
t = regexpi(tmp, 'K=\D*([0-9]+)', 'tokens');
if (length(t) > 0) size3 = str2num(t{1}{1}); else size3 = 0; end;       

%try
    while((length(tmp) >= 4) && (max(tmp ~= -1)) && strcmpi(tmp([1:4]), 'ZONE') && zone_number < n_zones)

        t = regexpi(tmp, 'T="(.*)"', 'tokens');      
        if (length(t) > 0) time = str2num(t{1}{1}); end;
        
%        t = regexpi(tmp, 'I=\D*([0-9]+)', 'tokens');
%        if (length(t) > 0) size1 = str2num(t{1}{1}); else size1 = 0; end;
%        t = regexpi(tmp, 'J=\D*([0-9]+)', 'tokens');
%        if (length(t) > 0) size2 = str2num(t{1}{1}); else size2 = 0; end;       
%        t = regexpi(tmp, 'K=\D*([0-9]+)', 'tokens');
%        if (length(t) > 0) size3 = str2num(t{1}{1}); else size3 = 0; end;       
        

%        if ((max(size(size1)) + max(size(size2)) + max(size(size2))) == 0)
        if ((size1 + size2 + size2) == 0)
            f_arr = fscanf(fid, '%f', [n_of_vars, inf]);
            size1 = size(f_arr, 2);
            size2 = 1;
            size3 = 1;
        else
%            if (max(size(size2)) == 0) size2 = 1; end;
%            if (max(size(size3)) == 0) size3 = 1; end;        
            if (size2 == 0) size2 = 1; end;
            if (size3 == 0) size3 = 1; end;        
            total_size = size1 * size2 * size3;
            f_arr = fscanf(fid, '%f', [n_of_vars, total_size]);
            if (size2 == 1 && size3 == 1) size1 = max(size(f_arr)); end;
        end

        tmp = fgetl(fid); % read to the end of line
        while (~feof(fid) && ((length(tmp) < 4) || (strcmpi(tmp([1:4]), 'ZONE') ~= 1)))
            tmp = fgetl(fid);
        end;

        if (zn >= first_zone-1) 
            if  sz <= 0
                if (zone_number > 0) 
                    fprintf('%s\n',['Loading: time = ', num2str(time), '...'] );
                    system('tput cuu 1');
                end;
                sz = skip_zones+1;
                zone_number = zone_number + 1;
                for var_number = 1:n_of_vars
                    if (length(time) > 0) vars(var_number).time(zone_number) = time; end; %stupid, but do not know other way to check if it is number...
                    vars(var_number).size1 = size1;
                    vars(var_number).size2 = size2;
                    vars(var_number).size3 = size3;
                    vars(var_number).arr(zone_number, :, :, :) = reshape(f_arr(var_number, :), [size1, size2, size3]); 
    %                vars(var_number).arr(zone_number, :, :, :) = reshape(f_arr(var_number, :), [size3, size2, size1]); 
                end
            end
            sz = sz - 1;
        else 
            zn = zn + 1;            
        end
        
    end
%catch
    
%end
    
    
if (zone_number == 0)

    f_arr = [sscanf(tmp, '%f', [n_of_vars, 1])'; fscanf(fid, '%f', [n_of_vars, inf])']'
    
    size1 = max(size(f_arr));
    size2 = 1;
    size3 = 1;

    for var_number = 1:n_of_vars
        vars(var_number).size1 = size1;
        vars(var_number).size2 = size2;
        vars(var_number).size3 = size3;
        vars(var_number).arr(:, :, :) = reshape(f_arr(var_number, :), [size1, size2, size3]); 
    end

end

if (dopermute) 
    for it = 1:n_of_vars
        vars(it).arr = permute(vars(it).arr, [1, 4, 3, 2]);
        var_tmp = vars(it).size1;
        vars(it).size1 = vars(it).size3;
        vars(it).size3 = var_tmp;
    end
end

if (dosqueeze) 
    for it = 1:n_of_vars
        vars(it).arr = squeeze(vars(it).arr);
    end    
end

nout = max(nargout,1);
for k=1:nout
    try 
        varargout(k) = {vars(k)};
    catch 
        err = lasterror;
        if (err.identifier == 'MATLAB:badsubscript')
            err.message = ['Only ', num2str(k-1), ' variables in the file (', num2str(nout), ' requested).'];
            rethrow(err);
        else 
            rethrow(lasterror)
        end
    end
end

fclose (fid);

return ;
