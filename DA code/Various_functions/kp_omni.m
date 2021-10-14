clear

date_start = datenum(2012,10,01);
date_end   = datenum(2013,10,01)+5;

[time, kp] = omniwebdata(date_start, date_end, 38);

T = time - time(1);
K = kp.arr ./ 10;

filename = './Kp_01Oct12_01Oct13.dat';

dlmwrite(filename,[T, K],'delimiter','\t');
