function [time, varargout] = omniwebdata(start_date, end_date, varargin)
    %% 
    % OMNIWEBDATA gets data from OMNIWEB (http://omniweb.gsfc.nasa.gov)
    %
    % [time, varargout] = omniwebdata(start_date, end_date, varargin)
    % 
    % Input:
    %   start_date and end_date in Matlab date format
    %   data arguments, vargin, as number, indicating what you whant to get
    %   (see below)
    %
    % Output:
    %   time - time
    %   varargout - data structures as:
    %     data.time - time for this data output
    %     data.arr - data values array
    % 
    % Example:
    %     omniwebdata(now-10, now, 4, 9, 38, 40);
    %
    % vargin values:
    %    3 - Bartels Rotation Number
    %    4 - IMF Spacecraft ID 
    %    5 - Plasma Spacecraft ID 
    %    6 - # Fine Scale Points in IMF Avgs 
    %    7 - # Fine Scale Points in Plasma Avgs
    %    8 - IMF Magnitude Avg, nT 
    %    9 - Magnitude, Avg IMF Vr, nT 
    %    10 - Lat. of Avg. IMF, deg. 
    %    11 - Long. of Avg. IMF, deg.  
    %    12 - Bx, GSE/GSM, nT 
    %    13 - By, GSE, nT  
    %    14 - Bz, GSE, nT 
    %    15 - By, GSM, nT 
    %    16 - Bz, GSM, nT 
    %    17 - Sigma in IMF Magnitude Avg. 
    %    18 - Sigma in IMF Vector Avg 
    %    19 - Sigma Bx, nT 
    %    20 - Sigma By, nT 
    %    21 - Sigma Bz, nT             
    %    22 - Proton Temperature, K 
    %    23 - Proton Density, n/cc 
    %    24 - Flow Speed, km/sec 
    %    25 - Flow Longitude, deg. 
    %    26 - Flow Latitude, deg. 
    %    27 - Alpha/Proton Density Ratio 
    %    28 - Flow Pressure, nPa 
    %    29 - Sigma-T 
    %    30 - Sigma-Np 
    %    31 - Sigma-V 
    %    32 - Sigma-Flow-Longitude 
    %    33 - Sigma-Flow-Latitude 
    %    34 - Sigma-Alpha/Proton Ratio 
    %    35 - Ey - Electric Field, mV/m 
    %    36 - Plasma Beta 
    %    37 - Alfven Mach Number 
    %    38 - Kp*10 Index 
    %    39 - R Sunspot Number 
    %    40 - Dst Index, nT 
    %    41 - AE Index, nT
    %    42 - Proton Flux* > 1 MeV 
    %    43 - Proton Flux* > 2 MeV 
    %    44 - Proton Flux* > 4 MeV     
    %    45 - Proton Flux* > 10 MeV 
    %    46 - Proton Flux* > 30 MeV 
    %    47 - Proton Flux* > 60 MeV 
    %    48 - Magnetospheric Flux Flag 
    %    49 - ap index, nT 
    %    50 - Solar index F10.7 
    %    51 - Polar Cap (PC) index from Thule     
    %    52 - AL Index, nT 
    %    53 - AU Index, nT 
    
    % Author: Dmitriy Subbotin
    % Email: subbotin@ucla.edu
    % Last change: 20014-01-01 by Dmitriy Subbotin
    % The work was supported by:
    %   LANL grant 12-LR-235337, PI Yuri Shprits
    %   NASA grant NNX09AF51G, PI Yuri Shprits
    
    vars = '';
    varnum = nargin - 2;
    if (varnum > 0)
        for it = 1:1:varnum
            vars = [vars, 'vars=', num2str(varargin{it}), '&'];
        end
    else
        'Error - no values requested'
        nargin
        varargin
        return
    end

 %
    PATH = getenv('PATH');
    setenv('PATH', [PATH ':/usr/local/bin']);
 %
    
    start_date_url = datestr(start_date, 'yyyymmdd');
    end_date_url = datestr(end_date, 'yyyymmdd');
    [~, ~] = system(['wget "http://omniweb.gsfc.nasa.gov/cgi/nx1.cgi?activity=retrieve&res=hour&spacecraft=omni2&start_date=', start_date_url, '&end_date=', end_date_url, '&maxdays=31&', vars, 'scale=Linear&view=0&nsum=1&paper=0&charsize=&xstyle=0&ystyle=0&symbol=0&symsize=&linestyle=solid&table=0&imagex=640&imagey=480&color=&back=" -O test_wget.txt']);
    fid = fopen('test_wget.txt');
    input_format = ['%d %d %d ', repmat('%f ', [1 varnum])];
    fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid); fgetl(fid); 
    for i = 1:varnum
        data(i).name = fgetl(fid); 
    end
    fgetl(fid); fgetl(fid);
    result = textscan(fid, input_format);
%    result = textscan(fid, input_format, 'headerLines', 5 + varnum + 2);
    fclose(fid);

    t_length = length(result{1});
    z_vec = zeros(t_length,1);
    year = double(result{1});
    doy = double(result{2});
    hh = double(result{3});
    
    time = datenum(year, z_vec, doy, hh, z_vec, z_vec);
    
    start_idx = min(find(time >= start_date));
    end_idx = min(find(time >= end_date));
    if isempty(end_idx)
        end_idx = length(time);
    end
    result_idx = [start_idx:end_idx];    
   
    time = time(result_idx);
    
    for i = 1:varnum
        data(i).arr = result{3+i};
        data(i).arr = data(i).arr(result_idx);
        data(i).time = time;
        varargout(i) = {data(i)};
    end
    
return
