function dim = find_DIM(yr,mo)
   dims=[31,28,31,30,31,30,31,31,30,31,30,31];
   if rem(yr,4) == 0 
       dims(2)=29;
   end
   dim=dims(mo);
return;
