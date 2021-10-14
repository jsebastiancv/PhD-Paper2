function [tstr] = format_time(t1);

% outputs a string of the current time elapsed
% input is output from toc command

ms = round(mod(t1,1)*100);
if t1 < 60
    sec = floor(t1);
    tstr = sprintf('00:00:%02i:%02i',round(t1),ms);
else
    sec = floor(mod(t1,60));
    mts = floor(mod(t1/60,60));
    hrs = floor(mod(t1/3600,24));
    if t1 < 86400
        tstr = sprintf('%02i:%02i:%02i:%02i',hrs,mts,sec,ms);
    else
        days = floor(t1/86400);
        tstr = sprintf('%02i dy, %02i:%02i:%02i:%02i',days,hrs,mts,sec,ms);
    end
end 

