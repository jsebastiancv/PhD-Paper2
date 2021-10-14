function [tstr] = estimation_time(start_tic,it,nt);

el_time = toc(start_tic); %seconds elapsed
it_time = el_time/it; %average iteration time
c_time = it_time*(nt-it); %estimated completion time

if c_time < 60
    sec = floor(c_time);
    tstr = sprintf('00:00:%02i',round(c_time));
else
    sec = floor(mod(c_time,60));
    mts = floor(mod(c_time/60,60));
    hrs = floor(mod(c_time/3600,24));
    if c_time < 86400
        tstr = sprintf('%02i:%02i:%02i',hrs,mts,sec);
    else
        days = floor(ctime/86400);
        tstr = sprintf('%02i dy, %02i:%02i:%02i',days,hrs,mts,sec);
    end
end 

