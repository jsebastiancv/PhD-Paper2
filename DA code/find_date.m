% kp_all = load('Kp_01Oct12_01Oct13.dat');
start_date = datenum('01-Apr-2013');
end_date = datenum('01-Aug-2013');
date_vec = linspace(start_date,end_date,122*24+1);

kp_all(:,1) = date_vec;
start_date_storm = datenum('15-Jun-2013,00:00');
% end_date_storm = datenum('15-Oct-2012');

idx_start_date_storm = find(date_vec == start_date_storm)
idx_end_date_storm = find(date_vec == end_date_storm);

kp_all([idx_start_date_storm:idx_end_date_storm],2)
