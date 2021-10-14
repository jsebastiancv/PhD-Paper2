function loop_counter(inx,nt,varargin);

use_tput=false;
modval=1;
pr_text=' ';
if nargin > 2
    for it=1:1:nargin-2
        if strcmp(varargin(it), 'tput')
            use_tput=true;
        elseif strcmp(varargin(it),'mod')
            modval = varargin{it+1};
            it = it + 1;
        elseif strcmp(varargin(it),'pretext')
            pr_text = varargin{it+1};
            it = it + 1;
        end
    end
end

if use_tput
    if mod(inx,modval) == 0
        system('tput cuu 1;tput el');
        fprintf('%s%s/%s\n',pr_text,num2str(inx),num2str(nt));
    end
else
    if mod(inx,modval) == 0
        fprintf('%s%s/%s\n',pr_text,num2str(inx),num2str(nt));
    end
end

