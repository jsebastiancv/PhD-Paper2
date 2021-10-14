% clc
% clear 
% close all
addpath([getenv('HOME'),'/Documents/VERB/22. 3D VERB DA_Spectrum_LT/Code/Various_functions/']);
addpath([getenv('HOME'),'/Documents/datalibrary/core/read/']);

%% 1. Load GOES data 

for iter = 1 : 3
    
    if iter == 1 
        sDate = datenum('01-Oct-2015');
        eDate = datenum('01-May-2016');
%         eDate = datenum('05-Nov-2012');

    elseif iter == 2
        sDate = datenum('01-May-2016');
        eDate = datenum('01-Jun-2016');
    elseif iter == 3
        sDate = datenum('01-Jun-2016');
        eDate = datenum('01-Aug-2016');
    end
    
    SatStruct.mission = 'goes'; 
    SatStruct.satellite = 'goes13'; 
    SatStruct.instrument = 'magedandepead';
    SatStruct.path = '/export/.automounter/mag5/data/rbm-data';

    SatStruct.mfm = 'TS07Dmid15';
    SatStruct.version = 'ver4';

    [Alphalocal13, Energy13, AlphaEqModel13, AlphaEqReal13, InvMu13, InvK13, InvMuReal13, ...
           Lstar13, Flux13, PSD13, Position13, MLT13, BVecModel13, BTotalModel13, Time13] = ...
                read_pmf(sDate, eDate, SatStruct,  'AlphaLocal', ...
                                                   'Energy', ...
                                                   'AlphaEqModel', ...
                                                   'AlphaEqReal', ...
                                                   'InvMu', ...
                                                   'InvK', ...
                                                   'InvMuReal', ...
                                                   'Lstar', ...
                                                   'Flux', ...
                                                   'PSD', ...
                                                   'Position', ...
                                                   'MLT' , ...
                                                   'BVecModel', ...
                                                   'BTotalModel', ...
                                                   'Time');

    SatStruct.satellite = 'goes15'; 

    [Alphalocal15, Energy15, AlphaEqModel15, AlphaEqReal15, InvMu15, InvK15, InvMuReal15, ...
           Lstar15, Flux15, PSD15, Position15, MLT15, BVecModel15, BTotalModel15, Time15] = ...
                read_pmf(sDate, eDate, SatStruct,  'AlphaLocal', ...
                                                   'Energy', ...
                                                   'AlphaEqModel', ...
                                                   'AlphaEqReal', ...
                                                   'InvMu', ...
                                                   'InvK', ...
                                                   'InvMuReal', ...
                                                   'Lstar', ...
                                                   'Flux', ...
                                                   'PSD', ...
                                                   'Position', ...
                                                   'MLT' , ...
                                                   'BVecModel', ...
                                                   'BTotalModel', ...
                                                   'Time');                                                            
    
    
    if iter == 1
        Flux13_A = Flux13; AlphaEqModel13_A = AlphaEqModel13; Lstar13_A = Lstar13; InvK13_A = InvK13; InvMu13_A = InvMuReal13; Time13_A = Time13;         
        Flux15_A = Flux15; AlphaEqModel15_A = AlphaEqModel15; Lstar15_A = Lstar15; InvK15_A = InvK15; InvMu15_A = InvMuReal15; Time15_A = Time15;                  
    elseif iter == 2
        Flux13_B = Flux15; AlphaEqModel13_B = AlphaEqModel15; Lstar13_B = Lstar15; InvK13_B = InvK15; InvMu13_B = InvMuReal15; Time13_B = Time15;        
%         Flux13_B = Flux13; AlphaEqModel13_B = AlphaEqModel13; Lstar13_B = Lstar13; InvK13_B = InvK13; InvMu13_B = InvMuReal13; Time13_B = Time13;        
        Flux15_B = Flux15; AlphaEqModel15_B = AlphaEqModel15; Lstar15_B = Lstar15; InvK15_B = InvK15; InvMu15_B = InvMuReal15; Time15_B = Time15;                  
    elseif iter == 3
        Flux13_C = Flux13; AlphaEqModel13_C = AlphaEqModel13; Lstar13_C = Lstar13; InvK13_C = InvK13; InvMu13_C = InvMuReal13; Time13_C = Time13;        
        Flux15_C = Flux15; AlphaEqModel15_C = AlphaEqModel15; Lstar15_C = Lstar15; InvK15_C = InvK15; InvMu15_C = InvMuReal15; Time15_C = Time15;                  
    end

    clearvars Flux13 AlphaEqModel13 Lstar13 InvK13 InvMu13 Time13 Flux15 AlphaEqModel15 Lstar15 InvK15 InvMu15 Time15               

end            

Flux13 = vertcat(Flux13_A.arr,Flux13_B.arr,Flux13_C.arr);
Flux15 = vertcat(Flux15_A.arr,Flux15_B.arr,Flux15_C.arr);

AlphaEqModel13 = vertcat(AlphaEqModel13_A.arr,AlphaEqModel13_B.arr,AlphaEqModel13_C.arr);
AlphaEqModel15 = vertcat(AlphaEqModel15_A.arr,AlphaEqModel15_B.arr,AlphaEqModel15_C.arr);

Lstar13 = vertcat(Lstar13_A.arr,Lstar13_B.arr,Lstar13_C.arr);
Lstar15 = vertcat(Lstar15_A.arr,Lstar15_B.arr,Lstar15_C.arr);

InvK13 = vertcat(InvK13_A.arr,InvK13_B.arr,InvK13_C.arr);
InvK15 = vertcat(InvK15_A.arr,InvK15_B.arr,InvK15_C.arr);

InvMu13 = vertcat(InvMu13_A.arr,InvMu13_B.arr,InvMu13_C.arr);
InvMu15 = vertcat(InvMu15_A.arr,InvMu15_B.arr,InvMu15_C.arr);

Time13 = vertcat(Time13_A.arr,Time13_B.arr,Time13_C.arr);
Time15 = vertcat(Time15_A.arr,Time15_B.arr,Time15_C.arr);

%% 2. Average Flux, AlphaEqModel and Lstar in intervals centered around the half hour

t_end = size(Time13);
% tvec_init = [1, linspace(7,8911,24*31-1), 8923]; % for 1 month
% tvec_end = [6, linspace(18,8922,24*31-1), 8928]; % for 1 month

tvec_init = [1, linspace(7,t_end(1)-17,24*(t_end(1)/(12*24))-1), t_end(1)-5];
tvec_end = [6, linspace(18,t_end(1)-6,24*(t_end(1)/(12*24))-1), t_end(1)];
s = size(tvec_end);

for t = 1:s(2)
    Flux13_bin(t,:,:) = nanmean(Flux13(tvec_init(t):tvec_end(t),:,:));
    AlphaEqModel13_bin(t,:) = nanmean(AlphaEqModel13(tvec_init(t):tvec_end(t),:));
    Lstar13_bin(t,:) = nanmean(Lstar13(tvec_init(t):tvec_end(t),:));
    InvK13_bin(t,:) = nanmean(InvK13(tvec_init(t):tvec_end(t),:));
    InvMu13_bin(t,:,:) = nanmean(InvMu13(tvec_init(t):tvec_end(t),:,:));
    
    
    Flux15_bin(t,:,:) = nanmean(Flux15(tvec_init(t):tvec_end(t),:,:));
    AlphaEqModel15_bin(t,:) = nanmean(AlphaEqModel15(tvec_init(t):tvec_end(t),:));
    Lstar15_bin(t,:) = nanmean(Lstar15(tvec_init(t):tvec_end(t),:));
    InvK15_bin(t,:) = nanmean(InvK15(tvec_init(t):tvec_end(t),:));
    InvMu15_bin(t,:,:) = nanmean(InvMu15(tvec_init(t):tvec_end(t),:,:));
end

FluxGOES_bin = nanmean(cat(4,Flux13_bin,Flux15_bin),4);
InvMuGOES_bin = nanmean(cat(4,InvMu13_bin,InvMu15_bin),4);
LstarGOES_bin = nanmean(cat(3,Lstar13_bin,Lstar15_bin),3);
AlphaEqModelGOES_bin = nanmean(cat(3,AlphaEqModel13_bin,AlphaEqModel15_bin),3);
InvKGOES_bin = nanmean(cat(3,InvK13_bin,InvK15_bin),3);

FluxGOES_bin = fillmissing(FluxGOES_bin,'linear');
InvMuGOES_bin = fillmissing(InvMuGOES_bin,'linear');
LstarGOES_bin = fillmissing(LstarGOES_bin,'linear');
AlphaEqModelGOES_bin = fillmissing(AlphaEqModelGOES_bin,'linear');
InvKGOES_bin = fillmissing(InvKGOES_bin,'linear');

%% 3. Load 3D VERB grid

gridFileName = ['perp_grid_Lma.plt'];
fprintf('Load grid from %s.\n', gridFileName);    
[L, epc, alpha, pc] = load_plt(gridFileName, 'squeeze', 'permute');

%% 4. Adiabatic transport 

for t = 1 : s(2)
    t
    for e = 1:17 % 17 energy channels in GOES grid

        e_goes = Energy13.arr(e);

        for a = 1:18 % 18 pitch angles in GOES grid

            a_goes = AlphaEqModelGOES_bin(t,a);
            Lstar_goes = LstarGOES_bin(t,a);

            mu = InvMuGOES_bin(t,e,a); % mu at original location
            K = InvKGOES_bin(t,a); % K at original location
% 
%             mu = pc2mu(Lstar_goes, pfunc(e_goes) , a_goes);
%             K = Lalpha2K(Lstar_goes, a_goes);
%         
            J = K2Jc(K,mu);

            alpha_ad(e,a) = LmuJc2alpha(6.6,mu,J); % alpha at L=6.6 conserving mu and K
            energy_ad = Lmualpha2pc(6.6,mu,alpha_ad(e,a));
            energy_mev_ad(e,a) = Kfunc(energy_ad); % energy at L = 6.6 conserving mu and K

        end
    end

    %  5. Interpolate from GOES grid to VERB grid

    E_VERB = squeeze(epc.arr(end,:,:)); 
    E_VERB_log = log10(squeeze(epc.arr(end,:,:)));
    A_VERB = squeeze(alpha.arr(end,:,:));

    E_GOES = squeeze(energy_mev_ad);
    A_GOES = squeeze(alpha_ad);

    Flux_t = squeeze(FluxGOES_bin(t,:,:));

    E_GOES2 = log10(E_GOES(:));
    A_GOES2 = A_GOES(:);
    Fluxt_2 = log10(Flux_t(:));
    F = scatteredInterpolant(E_GOES2,A_GOES2,Fluxt_2,'linear','none');

    Fq = 10.^(F(E_VERB_log,A_VERB));

    % 6. Find location of NaNs

    Fq_nan = isnan(Fq);
    Fq_nan_tot = sum(Fq_nan);
    C_first = find(Fq_nan_tot~=101, 1 ); % first column where there are any available fluxes, start extrapolating there

    % 7. Use J_L7 to extrapolate for E 

    x_jl7 = log10([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 2e-1, 1e0, 3e0, 1e1, 2e1]);
    y_jl7 = [3e9,  2e7,  9e6,  3e5,  1e4,  2e3, 7e0, 3e-3, 1e-3, 1e-3];

    for x = C_first : 91

        try
            C_column = isnan(Fq(:,x));
            R_first = find(C_column~=1, 1);
            R_last = find(C_column~=1, 1, 'last');

            for y = 1 : 101

                if isnan(Fq(y,x)) && y<R_first
                    JL7_targetE = interp1(x_jl7,y_jl7,E_VERB_log(y,x));
                    JL7_refE = interp1(x_jl7,y_jl7,E_VERB_log(R_first,x));
                    Fq(y,x) = Fq(R_first,x) * JL7_targetE / JL7_refE;
                elseif isnan(Fq(y,x)) && y>R_last
                    JL7_targetE = interp1(x_jl7,y_jl7,E_VERB_log(y,x));
                    JL7_refE = interp1(x_jl7,y_jl7,E_VERB_log(R_last,x));
                    Fq(y,x) = Fq(R_last,x) * JL7_targetE / JL7_refE;
                end

            end

        catch

        end

    end

    % 8. Use sine dependence to extrapolate for pitch angle

    for x = 1 : 101

        try
            R_row = isnan(Fq(x,:));
            C_first = find(R_row~=1, 1);
            C_last = find(R_row~=1, 1, 'last');

            for y = 1 : 91
                if isnan(Fq(x,y)) && y<C_first
                    Fq(x,y) = Fq(x,C_first) * sin(A_VERB(x,y)) / sin(A_VERB(x,C_first));
                elseif isnan(Fq(x,y)) && y>C_last
                    Fq(x,y) = Fq(x,C_last) * sin(A_VERB(x,y)) / sin(A_VERB(x,C_last));
                end
            end

        catch 

        end    

    end

    % 9. Convert flux to PSD

    PSD_BC = Fq ./ (pfunc(E_VERB).^2);

    % 10. Save PSD in file with grid

    PSD_Lupper(t,:,:) = PSD_BC;
    
    clear PSD_BC

end

