function [f,matmodel]=matsolve(A,B,C,f,varargin)
f=squeeze(f);
B=squeeze(B);
C=squeeze(C);
A=squeeze(A);

sol_meth = 'inverse'; % old standard
man_args = 4;
if (nargin > man_args)
    for it=1:1:nargin-man_args
        if strcmp(varargin(it),'solution_method')
            sol_meth = varargin{it+1};
            it = it+1;
        end
    end 
end

if size(f,2) >1 
    f=f';
end

% In VERB we solve the Fokker-Planck equation for PSD with 3 matrices, A, B, and C
% the method is as follows:
%   A*f[t+1] = B*f[t] + C
if strcmp(sol_meth,'inverse');
    % This is the old way uses inverse of A, this can be inaccurate...
    %=> f[t+1] = A'(B*f[t]) + A' + C
    iA=inv(squeeze(A));
    matmodel=iA*diag(B);
    f=matmodel*f+iA*C;
elseif strcmp(sol_meth,'tridiag')
    % set it up the same way as cpp code
%    RHS = diag(B)*f + C;
    RHS = B.*f + C; %elementwise multiplication and addition (no need to create a diagonal array)
    a = diag(A);
    c = a; %create same size
    b = a; % " "
    c(1:end-1) = diag(A,1); % upper diagonal
    b(2:end) = diag(A,-1);  % lower diagonal
    f = tridiag(a,b,c,RHS); % tridiag solution Thomas
    matmodel=tridiag(a,b,c,B);
elseif strcmp(sol_meth,'Gaussian')
%    RHS = diag(B)*f + C;
    RHS = B.*f + C; %elementwise multiplication and addition (no need to create a diagonal array)
    f = A\RHS;  
    matmodel=A\diag(B);
end 
