function [tstr] = now_time(varargin);

% default prints now, but can take a datenum input
time = now;

if nargin > 1
    error('only one argument is allowed');
elseif nargin == 1
    time = varargin{1};
end
% outputs a string of the specified time
tstr = datestr(time,'mmm-dd, HH:MM:SS');

