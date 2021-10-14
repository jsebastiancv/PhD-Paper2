
addpath([getenv('HOME'),'/Documents/datalibrary/core/read/']);

clc
clear
close all

sDate = datenum('01-Oct-2012');
eDate = datenum('01-Apr-2013');

%% Load satellite data

SatStruct.mission = 'goes'; 
SatStruct.satellite = 'goes13'; 
SatStruct.instrument = 'magedandepead';
SatStruct.path = '/export/.automounter/mag5/data/rbm-data';

SatStruct.mfm = 'TS07Dmid15';
SatStruct.version = 'ver4';

[Alphalocal, Energy, AlphaEqModel, AlphaEqReal, InvMu, InvK, InvMuReal, ...
                Lstar, Flux, PSD, Position, MLT, BVecModel, BTotalModel, Time] = ...
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

    
%% Interpolate to a target value of Mu and K

target_mu = 300;
target_k = 0.11;
psd_interp = interpolate_psd(log10(PSD.arr), InvMuReal.arr, InvK.arr, target_mu, target_k);

%% Plot interpolated data

subplot(3,1,1)
scatter(Time.arr, squeeze(Lstar.arr(:,end)), 8, psd_interp, 'filled');
ylabel(strcat('L*'))

axis tight
% ylim([3 6.6])
datetick('x')
title(strcat('\mu = ', num2str(target_mu),' MeV/G, K = ', num2str(target_k), ' G^{1/2} R_E'))

colorbar
colormap jet
ax = gca;

% ax.CLim = [-8 -4];

c1 = colorbar;
c1.Label.String = 'log_{10} PSD [(c/cm/MeV)^{3}]';
hold on
