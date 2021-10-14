function arr=load_matrix(filename)

%loads A/B/C matrix from .dat textfiles produced by VERB_MM_new
%combines A diagonals, returns 4d matrix for A
%returns 3d matrix for B & C

fid = fopen(filename, 'rt');

f_names = fgetl(fid);
while (isempty(f_names) || f_names(1) == '#')
    f_names = fgetl(fid);
end

f_names = regexpi(f_names, '"(.*?)"', 'match');
for it = 1:size(f_names,2)
    f_names{it} = strrep(f_names{it}, '"', ''); %reads in diagonal names for A matrix
end

tmp = fgetl(fid);
while length(tmp) < 2 || strcmpi(tmp(1:2), 'I=') ~= 1
    tmp = fgetl(fid);
end

t = regexpi(tmp, 'I=\D*([0-9]+)x([0-9]+)', 'tokens');
if ~isempty(t) 
    sizecalc2 = str2double(t{1}{2}); else sizecalc2 = 1; end; 
%sizecalc2 = 1 for B/C Matrix; sizecalc2 = # of diagonals for A matrix
t = regexpi(tmp, 'I=\D*([0-9]+)', 'tokens');
if (~isempty(t))
    sizecalc1 = str2double(t{1}{1}); else sizecalc1 = 1; end;
t = regexpi(tmp, 'J=\D*([0-9]+)', 'tokens');
if (~isempty(t))
    size2 = str2double(t{1}{1}); else size2 = 1; end;       
t = regexpi(tmp, 'K=\D*([0-9]+)', 'tokens');
if (~isempty(t))
    size3 = str2double(t{1}{1}); else size3 = 1; end;       

if size(f_names, 2) ~= sizecalc2
    error('ERROR: NUMBER OF DIAGONALS INCORRECT IN CALCULATION MATRIX');
end
if sizecalc2 == 1 % B/C Matrix
    lastdim = 1;
else
    lastdim = sizecalc1; % A Matrix
end
arr = nan(size2, size3, sizecalc1, lastdim);
i = 1; %calcmatrix iterator (sizecalc2)
j = 1; %x of calcmatrix4d (size2)
k = 1; %y of calcmatrix4d (size3)
Adiags = zeros(lastdim, lastdim);
while ~feof(fid)
    A = [];
    while isempty(A) && ~feof(fid)
        A = fscanf(fid, '%f', [1 sizecalc1]); %reads line into A
        fgetl(fid);
    end
    
    if feof(fid) && isempty(A)
        fclose(fid);
        return
    end
    if sizecalc2 > 1  
        d_n = str2double(f_names{i}); %gets correct diagonal from first line
        if (d_n <= 0)
            dg = diag(A(1-d_n:end), d_n);
        else
            dg = diag(A(1:end-d_n), d_n);
        end
        
        Adiags = Adiags + dg; %converts to diagonal
        i = i+1;
        if i == sizecalc2+1
            i = 1;
            arr(j, k, :, :) = Adiags;
            k = k+1;
            Adiags = zeros(sizecalc1, sizecalc1);
        end
    else
        arr(j, k, :) = A;
        k = k+1;
    end
    if mod(k-1, size3) == 0 && k ~= 1
        j = j+1;
        k = 1;
    end
end

fclose(fid);
