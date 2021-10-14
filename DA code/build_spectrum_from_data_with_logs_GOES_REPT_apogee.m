% clc
% close all
% clear

addpath([getenv('HOME'),'/Documents/VERB/3D VERB DA V4/Code/Various_functions/']);
addpath([getenv('HOME'),'/Documents/datalibrary/core/read/']);

%% 1.1 Load GOES data
sDate = datenum('01-Mar-2013');
eDate = datenum('15-Apr-2013');
SatStruct.mission = 'goes'; 
SatStruct.satellite = 'goes13'; 
SatStruct.instrument = 'magedandepead';
SatStruct.path = '/export/.automounter/mag5/rbm/data/rbm-data';

SatStruct.mfm = 'Ts07d';
SatStruct.version = 'ver3';

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

%% 1.2 Load rept data
sDate = datenum('01-Mar-2013');
eDate = datenum('15-Apr-2013');
SatStruct.mission = 'rbsp'; 
SatStruct.satellite = 'rbspa'; 
SatStruct.instrument = 'rept';
SatStruct.path = '/export/.automounter/mag5/rbm/data/rbm-data';

SatStruct.mfm = 'Ts07d';
SatStruct.version = 'ver3';

[AlphalocalA, EnergyA, AlphaEqModelA, AlphaEqRealA, InvMuA, InvKA, InvMuRealA, ...
       LstarA, FluxA, PSDA, PositionA, MLTA, BVecModelA, BTotalModelA, TimeA] = ...
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
                                           
SatStruct.satellite = 'rbspb'; 

[AlphalocalB, EnergyB, AlphaEqModelB, AlphaEqRealB, InvMuB, InvKB, InvMuRealB, ...
       LstarB, FluxB, PSDB, PositionB, MLTB, BVecModelB, BTotalModelB, TimeB] = ...
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

%% 1.3 Compute R for RBSP A and B

R_RBSP_A = sqrt(PositionA.arr(:,1).^2 + PositionA.arr(:,2).^2 + PositionA.arr(:,3).^2);
R_RBSP_B = sqrt(PositionB.arr(:,1).^2 + PositionB.arr(:,2).^2 + PositionB.arr(:,3).^2);

%% 2. Average Flux, AlphaEqModel and Lstar in intervals centered around the half hour

t_end = size(Time13.arr);
% tvec_init = [1, linspace(7,8911,24*31-1), 8923]; % for 1 month
% tvec_end = [6, linspace(18,8922,24*31-1), 8928]; % for 1 month

tvec_init = [1, linspace(7,t_end(1)-17,24*(t_end(1)/(12*24))-1), t_end(1)-5];
tvec_end = [6, linspace(18,t_end(1)-6,24*(t_end(1)/(12*24))-1), t_end(1)];
s = size(tvec_end);

for t = 1:s(2)
    Flux13_bin(t,:,:) = nanmean(Flux13.arr(tvec_init(t):tvec_end(t),:,:));
    AlphaEqModel13_bin(t,:) = nanmean(AlphaEqModel13.arr(tvec_init(t):tvec_end(t),:));
    Lstar13_bin(t,:) = nanmean(Lstar13.arr(tvec_init(t):tvec_end(t),:));
    InvK13_bin(t,:) = nanmean(InvK13.arr(tvec_init(t):tvec_end(t),:));
    InvMu13_bin(t,:,:) = nanmean(InvMu13.arr(tvec_init(t):tvec_end(t),:,:));
       
    Flux15_bin(t,:,:) = nanmean(Flux15.arr(tvec_init(t):tvec_end(t),:,:));
    AlphaEqModel15_bin(t,:) = nanmean(AlphaEqModel15.arr(tvec_init(t):tvec_end(t),:));
    Lstar15_bin(t,:) = nanmean(Lstar15.arr(tvec_init(t):tvec_end(t),:));
    InvK15_bin(t,:) = nanmean(InvK15.arr(tvec_init(t):tvec_end(t),:));
    InvMu15_bin(t,:,:) = nanmean(InvMu15.arr(tvec_init(t):tvec_end(t),:,:));
    
    FluxA_bin(t,:,:) = nanmean(FluxA.arr(tvec_init(t):tvec_end(t),:,:));
    AlphaEqModelA_bin(t,:) = nanmean(AlphaEqModelA.arr(tvec_init(t):tvec_end(t),:));
    LstarA_bin(t,:) = nanmean(LstarA.arr(tvec_init(t):tvec_end(t),:));
    InvKA_bin(t,:) = nanmean(InvKA.arr(tvec_init(t):tvec_end(t),:));
    InvMuA_bin(t,:,:) = nanmean(InvMuA.arr(tvec_init(t):tvec_end(t),:,:));
    RA_bin(t) = nanmean(R_RBSP_A(tvec_init(t):tvec_end(t)));
    
    FluxB_bin(t,:,:) = nanmean(FluxB.arr(tvec_init(t):tvec_end(t),:,:));
    AlphaEqModelB_bin(t,:) = nanmean(AlphaEqModelB.arr(tvec_init(t):tvec_end(t),:));
    LstarB_bin(t,:) = nanmean(LstarB.arr(tvec_init(t):tvec_end(t),:));
    InvKB_bin(t,:) = nanmean(InvKB.arr(tvec_init(t):tvec_end(t),:));
    InvMuB_bin(t,:,:) = nanmean(InvMuB.arr(tvec_init(t):tvec_end(t),:,:));
    RB_bin(t) = nanmean(R_RBSP_B(tvec_init(t):tvec_end(t)));
    
end

%% 2.1 Delete (make NaNs) measurements with R < 5.5

for t = 1:s(2)

    if RA_bin(t) < 5.5
        FluxA_bin(t,:,:) = NaN;
        AlphaEqModelA_bin(t,:) = NaN;
        LstarA_bin(t,:) = NaN;
        InvKA_bin(t,:) = NaN;
        InvMuA_bin(t,:,:) = NaN;
    end

    if RB_bin(t) < 5.5
        FluxB_bin(t,:,:) = NaN;
        AlphaEqModelB_bin(t,:) = NaN;
        LstarB_bin(t,:) = NaN;
        InvKB_bin(t,:) = NaN;
        InvMuB_bin(t,:,:) = NaN;
    end
end


%%

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
%%
FluxRBSP_bin = nanmean(cat(4,FluxA_bin,FluxB_bin),4);
InvMuRBSP_bin = nanmean(cat(4,InvMuA_bin,InvMuB_bin),4);
LstarRBSP_bin = nanmean(cat(3,LstarA_bin,LstarB_bin),3);
AlphaEqModelRBSP_bin = nanmean(cat(3,AlphaEqModelA_bin,AlphaEqModelB_bin),3);
InvKRBSP_bin = nanmean(cat(3,InvKA_bin,InvKB_bin),3);

FluxRBSP_bin = fillmissing(FluxRBSP_bin,'linear');
InvMuRBSP_bin = fillmissing(InvMuRBSP_bin,'linear');
LstarRBSP_bin = fillmissing(LstarRBSP_bin,'linear');
AlphaEqModelRBSP_bin = fillmissing(AlphaEqModelRBSP_bin,'linear');
InvKRBSP_bin = fillmissing(InvKRBSP_bin,'linear');

%% 3. Load 3D VERB grid

gridFileName = ['perp_grid_Lma.plt'];
fprintf('Load grid from %s.\n', gridFileName);    
[L, epc, alpha, pc] = load_plt(gridFileName, 'squeeze', 'permute');

ech_goes = 9; % energy channel number 9 = 2 MeV
  
%% 4. Adiabatic transport of GOES measurements with energy <= 2 MeV

for t = 1 : s(2)  

    for e = 1:ech_goes 

        e_goes = Energy13.arr(e);

        for a = 1:18 % 18 pitch angles in GOES grid

            a_goes = AlphaEqModelGOES_bin(t,a);
            Lstar_goes = LstarGOES_bin(t,a);

        mu = InvMuGOES_bin(t,e,a); % mu at original location
        K = InvKGOES_bin(t,a); % K at original location

%         mu = pc2mu(Lstar_goes, pfunc(e_goes) , a_goes);
%         K = Lalpha2K(Lstar_goes, a_goes);
%            

            J = K2Jc(K,mu);

            alpha_ad_goes(e,a) = LmuJc2alpha(6.6,mu,J); % alpha at L=6.6 conserving mu and K
            energy_ad_goes = Lmualpha2pc(6.6,mu,alpha_ad_goes(e,a));
            energy_mev_ad_goes(e,a) = Kfunc(energy_ad_goes); % energy at L = 6.6 conserving mu and K

        end
    end


%  5. Interpolate from GOES grid to VERB grid up to 1.5 MeV (rows 1 to 67 in VERB grid)

    ech_verb = 67;

    E_VERB_low = squeeze(epc.arr(end,1:ech_verb,:)); 
    E_VERB_low_log = log10(squeeze(epc.arr(end,1:ech_verb,:)));
    A_VERB_low = squeeze(alpha.arr(end,1:ech_verb,:));

    E_GOES = squeeze(energy_mev_ad_goes);
    A_GOES = squeeze(alpha_ad_goes);

    Flux_t = squeeze(FluxGOES_bin(t,1:ech_goes,:));

    E_GOES2 = log10(E_GOES(:));
    A_GOES2 = A_GOES(:);
    Fluxt_2 = log10(Flux_t(:));
    F = scatteredInterpolant(E_GOES2,A_GOES2,Fluxt_2,'linear','none');

    Fq_low = 10.^(F(E_VERB_low_log,A_VERB_low));

%   6. Adiabatic transport of RBSP measurements with energy > 2 MeV

    ech_rbsp = 12; % 

    for e = 1:ech_rbsp 

        e_rbsp = EnergyA.arr(e);

        for a = 1:9 % 9 pitch angles in RBSP grid

            a_rbsp = AlphaEqModelRBSP_bin(t,a);
            Lstar_rbsp = LstarRBSP_bin(t,a);

        mu = InvMuRBSP_bin(t,e,a); % mu at original location
        K = InvKRBSP_bin(t,a); % K at original location
% 
%         mu = pc2mu(Lstar_rbsp, pfunc(e_rbsp) , a_rbsp);
%         K = Lalpha2K(Lstar_rbsp, a_rbsp);

            J = K2Jc(K,mu);

            alpha_ad_rbsp(e,a) = LmuJc2alpha(6.6,mu,J); % alpha at L=6.6 conserving mu and K
            energy_ad_rbsp = Lmualpha2pc(6.6,mu,alpha_ad_rbsp(e,a));
            energy_mev_ad_rbsp(e,a) = Kfunc(energy_ad_rbsp); % energy at L = 6.6 conserving mu and K

        end
    end    
    
%   7. Interpolate from RBSP grid to VERB grid from 1.5 MeV to 10 MeV (rows 68 to 101 in VERB grid)

    clear Flux_t Fluxt_2 F

    E_VERB_high = squeeze(epc.arr(end,ech_verb+1:end,:)); 
    E_VERB_high_log = log10(squeeze(epc.arr(end,ech_verb+1:end,:)));
    A_VERB_high = squeeze(alpha.arr(end,ech_verb+1:end,:));

    E_RBSP = squeeze(energy_mev_ad_rbsp);
    A_RBSP = squeeze(alpha_ad_rbsp);

    Flux_t = squeeze(FluxRBSP_bin(t,1:ech_rbsp,:));

    E_RBSP2 = log10(E_RBSP(:));
    A_RBSP2 = A_RBSP(:);
    Fluxt_2 = log10(Flux_t(:));
    F = scatteredInterpolant(E_RBSP2,A_RBSP2,Fluxt_2,'linear','none');

    Fq_high = 10.^(F(E_VERB_high_log,A_VERB_high));
    
%   8. Merge both interpolated fluxes
    Fq = [Fq_low;Fq_high];    
    
%   9. Find location of NaNs

    Fq_nan = isnan(Fq);
    Fq_nan_tot = sum(Fq_nan);
    C_first = find(Fq_nan_tot~=101, 1 ); % first column where there are any available fluxes, start extrapolating there

    E_VERB = squeeze(epc.arr(end,:,:)); 
    E_VERB_log = log10(squeeze(epc.arr(end,:,:)));
    A_VERB = squeeze(alpha.arr(end,:,:)); 

%    10. Use J_L7 to extrapolate for E 

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

%   11. Use sine dependence to extrapolate for pitch angle

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

%   12. Convert flux to PSD

    PSD_BC = Fq ./ (pfunc(E_VERB).^2);

%   13. Save PSD in file with grid

    PSD_Lupper(t,:,:) = PSD_BC;
    
    clear PSD_BC

end

