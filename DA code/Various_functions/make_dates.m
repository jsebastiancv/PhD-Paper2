function dates=make_dates(yr,mo)
   dim=find_DIM(yr,mo);
   dates=strcat(datestr(datenum(yr, mo,1),'yyyymmdd'),...
       'to',datestr(datenum(yr,mo,dim),'yyyymmdd'));
return

