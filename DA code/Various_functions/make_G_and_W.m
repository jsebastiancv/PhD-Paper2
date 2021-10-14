%%%% Written by Adam Kellerman

%%make_G_and_W, lsw_dat% B_X to AE
%format long;

%lsw_data format
% 1  -  year
% 2  -  doy
% 3  -  month
% 4  -  day
% 5  -  hour
% 6  -  minute
% 7  -  second
% 8  -  total time in second since Jan 1, 1970
% 9  -  altitude in km
% 10 - geographic latitude
% 11 - geographic longitude  (GDZ coordinate)
% 12 - BxIMF
% 13 - ByIMF
% 14 - BzIMF
% 15 - SWSpeed
% 16 - SWDensity
% 17 - SWPressure
% 18 - Dst_tmp
% 19 - Kp
% 20 - AE
% 21 - 23 G1 - G3
% 24 - 29 W1 - W6
clear G1
clear G2
clear G3
clear W1
clear W2
clear W3
clear W4
clear W5
clear W6

% Year=SW(:,1);
% Doy=SW(:,2);
% Hour=SW(:,3);
% Min=SW(:,4);

By=sw_data(:,18);
Bz=sw_data(:,19);
Vsw=sw_data(:,22);
Nsw=sw_data(:,26);
Pdyn=sw_data(:,28);
%AE=sw_data(:,38);
%Kp=lsw_dat(:,19);
%Dst_tmp=lsw_dat(:,18);
[ss1,tt1]=size(sw_data);
%%
% parameters for T2001 models
G1(1:ss1,1)=NaN;
G2(1:ss1,1)=NaN;
G3(1:ss1,1)=NaN;

% parameters for T2010 storm model
W1(1:ss1,1)=NaN;
W2(1:ss1,1)=NaN;
W3(1:ss1,1)=NaN;
W4(1:ss1,1)=NaN;
W5(1:ss1,1)=NaN;
W6(1:ss1,1)=NaN;

% for T2001 models
clear Bperp
clear h_bperp
clear ind1
clear theta
clear Bs
clear ind2
clear ind3

Bperp(:,1)=sqrt(Bz.*Bz+By.*By);
h_bperp(:,1)=(Bperp./40.0).*(Bperp./40.0)./(1.0+Bperp./40.0);
ind1=find(Bz>0);
theta(ind1,1)=atan(abs(By(ind1))./Bz(ind1));
Bs(ind1,1)=0.0;
ind2=find(Bz==0);
theta(ind2,1)=0.5;
Bs(ind2,1)=0.0;
ind3=find(Bz<0);
theta(ind3,1)=pi-atan(abs(By(ind3))./Bz(ind3));
Bs(ind3,1)=-Bz(ind3);
if ~isempty(theta)
    if length(theta) == length(Bz)
        
        
        
        G1_series=Vsw.*h_bperp.*sin(theta./2.0).*sin(theta./2.0).*sin(theta./2.0);
        G2_series=Vsw.*Bs./200.0;
        G3_series=Nsw.*Vsw.*Bs./2000.0;
                       
        % for T2010 storm model
        r1=0.39; d1=0.39; b1=0.80; g1=0.87;
        r2=0.70; d2=0.46; b2=0.18; g2=0.67;
        r3=0.031;d3=0.39; b3=2.32; g3=1.32;
        r4=0.58; d4=0.42; b4=1.25; g4=1.29;
        r5=1.15; d5=0.41; b5=1.60; g5=0.69;
        r6=0.88; d6=1.29; b6=2.40; g6=0.53;
        
        S1=(Nsw./5.0).^d1.*(Vsw./400.0).^b1.*(Bs./5.0).^g1;
        S2=(Nsw./5.0).^d2.*(Vsw./400.0).^b2.*(Bs./5.0).^g2;
        S3=(Nsw./5.0).^d3.*(Vsw./400.0).^b3.*(Bs./5.0).^g3;
        S4=(Nsw./5.0).^d4.*(Vsw./400.0).^b4.*(Bs./5.0).^g4;
        S5=(Nsw./5.0).^d5.*(Vsw./400.0).^b5.*(Bs./5.0).^g5;
        S6=(Nsw./5.0).^d6.*(Vsw./400.0).^b6.*(Bs./5.0).^g6;
        
        %%
        
        parfor i=61:ss1
            
            % calculate G1 - averaging over the preceding 1-hour interval
            % (Tsyganenko, 2010, JGR, 107(A8), 1176, equation (1))
            G1(i,1)=nanmean(G1_series((i-60):(i)));
            % calculate G2 - averaging over the preceding 1-hour interval
            % (Tsyganenko, 2010, JGR, 107(A8), 1176, equation (2) or Tsyganenko et al., 2010, JGR, 108(A5), 1209, below equation (1))
            G2(i,1)=nanmean(G2_series((i-60):(i)));
            % calculate G3 - averaging over the preceding 1-hour interval
            % (Tsyganenko et al., 2010, JGR, 108(A5), 1209, below equation (2))
            G3(i,1)=nanmean(G3_series((i-60):(i)));
                   
            
            % calculate W1 - averaging over the preceding 1-hour interval
            % (Tsyganenko and Sitnov, 2010, JGR, 110, A03208, equations (7) and (8))
            W1(i,1)=r1/12.0*(S1(i-55)*exp(r1/60.0*(-55.0))+S1(i-50)*exp(r1/60.0*(-50.0))+S1(i-45)*exp(r1/60.0*(-45.0))+...
                S1(i-40)*exp(r1/60.0*(-40.0))+S1(i-35)*exp(r1/60.0*(-35.0))+S1(i-30)*exp(r1/60.0*(-30.0))+...
                S1(i-25)*exp(r1/60.0*(-25.0))+S1(i-20)*exp(r1/60.0*(-20.0))+S1(i-15)*exp(r1/60.0*(-15.0))+...
                S1(i-10)*exp(r1/60.0*(-10.0))+S1(i-5)*exp(r1/60.0*(-5.0))+S1(i-0)*exp(r1/60.0*(-0.0)));
            % calculate W2 - averaging over the preceding 1-hour interval
            % (Tsyganenko and Sitnov, 2010, JGR, 110, A03208, equations (7) and (8))
            W2(i)=r2/12.0*(S2(i-55)*exp(r2/60.0*(-55.0))+S2(i-50)*exp(r2/60.0*(-50.0))+S2(i-45)*exp(r2/60.0*(-45.0))+...
                S2(i-40)*exp(r2/60.0*(-40.0))+S2(i-35)*exp(r2/60.0*(-35.0))+S2(i-30)*exp(r2/60.0*(-30.0))+...
                S2(i-25)*exp(r2/60.0*(-25.0))+S2(i-20)*exp(r2/60.0*(-20.0))+S2(i-15)*exp(r2/60.0*(-15.0))+...
                S2(i-10)*exp(r2/60.0*(-10.0))+S2(i-5)*exp(r2/60.0*(-5.0))+S2(i-0)*exp(r2/60.0*(-0.0)));
            % calculate W3 - averaging over the preceding 1-hour interval
            % (Tsyganenko and Sitnov, 2010, JGR, 110, A03208, equations (7) and (8))
            W3(i)=r3/12.0*(S3(i-55)*exp(r3/60.0*(-55.0))+S3(i-50)*exp(r3/60.0*(-50.0))+S3(i-45)*exp(r3/60.0*(-45.0))+...
                S3(i-40)*exp(r3/60.0*(-40.0))+S3(i-35)*exp(r3/60.0*(-35.0))+S3(i-30)*exp(r3/60.0*(-30.0))+...
                S3(i-25)*exp(r3/60.0*(-25.0))+S3(i-20)*exp(r3/60.0*(-20.0))+S3(i-15)*exp(r3/60.0*(-15.0))+...
                S3(i-10)*exp(r3/60.0*(-10.0))+S3(i-5)*exp(r3/60.0*(-5.0))+S3(i-0)*exp(r3/60.0*(-0.0)));
            % calculate W4 - averaging over the preceding 1-hour interval
            % (Tsyganenko and Sitnov, 2010, JGR, 110, A03208, equations (7) and (8))
            W4(i)=r4/12.0*(S4(i-55)*exp(r4/60.0*(-55.0))+S4(i-50)*exp(r4/60.0*(-50.0))+S4(i-45)*exp(r4/60.0*(-45.0))+...
                S4(i-40)*exp(r4/60.0*(-40.0))+S4(i-35)*exp(r4/60.0*(-35.0))+S4(i-30)*exp(r4/60.0*(-30.0))+...
                S4(i-25)*exp(r4/60.0*(-25.0))+S4(i-20)*exp(r4/60.0*(-20.0))+S4(i-15)*exp(r4/60.0*(-15.0))+...
                S4(i-10)*exp(r4/60.0*(-10.0))+S4(i-5)*exp(r4/60.0*(-5.0))+S4(i-0)*exp(r4/60.0*(-0.0)));
            % calculate W5 - averaging over the preceding 1-hour interval
            % (Tsyganenko and Sitnov, 2010, JGR, 110, A03208, equations (7) and (8))
            W5(i)=r5/12.0*(S5(i-55)*exp(r5/60.0*(-55.0))+S5(i-50)*exp(r5/60.0*(-50.0))+S5(i-45)*exp(r5/60.0*(-45.0))+...
                S5(i-40)*exp(r5/60.0*(-40.0))+S5(i-35)*exp(r5/60.0*(-35.0))+S5(i-30)*exp(r5/60.0*(-30.0))+...
                S5(i-25)*exp(r5/60.0*(-25.0))+S5(i-20)*exp(r5/60.0*(-20.0))+S5(i-15)*exp(r5/60.0*(-15.0))+...
                S5(i-10)*exp(r5/60.0*(-10.0))+S5(i-5)*exp(r5/60.0*(-5.0))+S5(i-0)*exp(r5/60.0*(-0.0)));
            % calculate W6 - averaging over the preceding 1-hour interval
            % (Tsyganenko and Sitnov, 2010, JGR, 110, A03208, equations (7) and (8))
            W6(i)=r6/12.0*(S6(i-55)*exp(r6/60.0*(-55.0))+S6(i-50)*exp(r6/60.0*(-50.0))+S6(i-45)*exp(r6/60.0*(-45.0))+...
                S6(i-40)*exp(r6/60.0*(-40.0))+S6(i-35)*exp(r6/60.0*(-35.0))+S6(i-30)*exp(r6/60.0*(-30.0))+...
                S6(i-25)*exp(r6/60.0*(-25.0))+S6(i-20)*exp(r6/60.0*(-20.0))+S6(i-15)*exp(r6/60.0*(-15.0))+...
                S6(i-10)*exp(r6/60.0*(-10.0))+S6(i-5)*exp(r6/60.0*(-5.0))+S6(i-0)*exp(r6/60.0*(-0.0)));
            
            if((mod(i,10000)==0) || i == ss1)
                fprintf('%s\n',[num2str(i),' of ', num2str(ss1)]);
            end
            
        end
        
%         lsw_dat(:,21)=G1;
%         lsw_dat(:,22)=G2;
%         lsw_dat(:,23)=G3;
%         lsw_dat(:,24)=W1;
%         lsw_dat(:,25)=W2;
%         lsw_dat(:,26)=W3;
%         lsw_dat(:,27)=W4;
%         lsw_dat(:,28)=W5;
%         lsw_dat(:,29)=W6;
    end
end
%%
