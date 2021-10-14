%
% Save arrays into .plt format
% save_plt(filename, arrays... )
%
%

function [varargout]=save_plt(filename, varargin)

%total size calculating
total_size = size(varargin{1}.arr, 1)*size(varargin{1}.arr, 2)*size(varargin{1}.arr, 3);

%adding all arrays to the 1 array
if (nargin > 1)
    %reshape first array to the 1d array and add to the total array
    total_array = reshape(permute(varargin{1}.arr, [3,2,1]), [total_size, 1]);
    n_of_vars = 1;
    % creating output format
    format_string = '%e\t';
    for it = 2:nargin-1
        %reshape others arrays to the 1d array and add to the total array
        total_array = [total_array, reshape(permute(varargin{it}.arr, [3,2,1]), [total_size, 1])];
        n_of_vars = n_of_vars + 1;
        format_string = [format_string, '%e\t'];
    end
end
format_string = [format_string, '\n'];

%remove NaN
total_array (find(~isfinite(total_array)))=1e-22; 

fid = fopen(filename, 'w');

%writing comments
for it = 1:nargin-1
    try
    for comment = varargin{it}.comments
        fprintf(fid, '%s\n', comment{1});
    end
    catch
    end
end
%    try
%        %trying to write comments, if exists
%        fprintf(fid, '%s\n', varargin{it}.comments);
%    catch
%        %else - nothing
%    end
%end

    
%writing header
fprintf(fid, 'VARIABLES = ');
for it = 1:nargin-1
    try
        %trying to write var-name, if exists
        fprintf(fid, '"%s", ', varargin{it}.name); % ...name{1}
    catch
        fprintf(fid, '"func.", ');
    end
end
fprintf(fid, '\nZONE T="..." I=%d, J=%d, K=%d\n', size(varargin{1}.arr, 3), size(varargin{1}.arr, 2), size(varargin{1}.arr, 1));

fprintf(fid, format_string, total_array');    

fclose(fid);
