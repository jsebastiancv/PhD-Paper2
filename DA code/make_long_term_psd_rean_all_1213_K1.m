clear 
clc
%close all

%% Settings
addpath([getenv('HOME'),'/Documents/VERB/22. 3D VERB DA_Spectrum_LT/Code/Various_functions/']);
addpath([getenv('HOME'),'/Documents/datalibrary/core/read/']);

mu_vec = 1000*[2,3,4,5,6,7,8,9,10,12,15];
    
for m = 1 : 11

    target_K = 0.1; % G^0.5 Re
    target_mu = mu_vec(m) % MeV/G
    angle = pi/180;
    mc2 = 0.511;
    fastint = true; % will find nearest values only
    psd_LT = [];
    Kp_LT = [];

    for i = 1% : 12
        clear sDate eDate strdate fileName data alpha energy K mu L PSD Kp sz PSD_new flux time
        target_K = 0.1; % G^0.5 Re
        target_mu = mu_vec(m); % MeV/G
        if i == 1; sDate = datenum('01-Oct-2012'); eDate = datenum('01-Nov-2012'); end
        if i == 2; sDate = datenum('01-Nov-2012'); eDate = datenum('01-Dec-2012'); end
        if i == 3; sDate = datenum('01-Dec-2012'); eDate = datenum('01-Jan-2013'); end
        if i == 4; sDate = datenum('01-Jan-2013'); eDate = datenum('01-Feb-2013'); end
        if i == 5; sDate = datenum('01-Feb-2013'); eDate = datenum('01-Mar-2013'); end
        if i == 6; sDate = datenum('01-Mar-2013'); eDate = datenum('01-Apr-2013'); end
        if i == 7; sDate = datenum('01-Apr-2013'); eDate = datenum('01-May-2013'); end
        if i == 8; sDate = datenum('01-May-2013'); eDate = datenum('01-Jun-2013'); end
        if i == 9; sDate = datenum('01-Jun-2013'); eDate = datenum('01-Jul-2013'); end
        if i == 10; sDate = datenum('01-Jul-2013'); eDate = datenum('01-Aug-2013'); end
        if i == 11; sDate = datenum('01-Aug-2013'); eDate = datenum('01-Sep-2013'); end
        if i == 12; sDate = datenum('01-Sep-2013'); eDate = datenum('01-Oct-2013'); end

        strdate = datestr(sDate, 'yyyymm');
        mfm = 'TS07Dmid15';

        % Load reanalysis file
%         fileName = ['reanalysis_TS07Dmid15_final/ReanalysisErr_EQE_LatestVERB_noMP_',strdate,'_Gaussian_onera_',mfm,'.mat'];
%         dataMod = load(fileName);
    %     
        fileName = ['reanalysis_TS07Dmid15_final_q=5000r_newDxx_RUN7/Reanalysis_EQE_LatestVERB_noMP_',strdate,'_Gaussian_onera_',mfm,'.mat'];
        data = load(fileName);
        alpha = data.SimInvAlpha;
        energy = data.SimInvEnergy;
        K = data.SimInvK;
        mu = data.SimInvMu;
        L = data.SimL;

        time = data.SimTime;
%         PSD = data.PSD_obs;
%         PSD = data.SimPSD;
        PSD = dataMod.ModPSD;

        Kp = data.Kp;
        sz = size(PSD);

        PSD_new = nan(sz(1),sz(2));

        % Interpolation
        if fastint
            fprintf('warning, using coarse interpolation\n')
            linx = floor(size(L,1)/2); % we need to do this on L,mu,alpha grid (approx at middle of grid)
            [~,kinx] = min(abs(target_K - K(linx,1,:)));
            tmpMu = squeeze(mu(:,:,kinx));
            [~,minx] = min(abs(target_mu - tmpMu(linx,:)));
            PSD_new = squeeze(PSD(:,:,minx,kinx));

            % Redefine values
            fprintf('redefining target mu and k to match nearest values\n')
            target_mu = round(tmpMu(linx,minx)); % mu is the same at all L
            target_K = K(linx,1,kinx); %K varies due to fixed alpha
        else
            fprintf('Interpolating reanalysis (this may take some time)...\n');
            t1=tic;

            PSD = log10(PSD);
            K = K.^(1/3);
            mu = log10(mu);

            for it=1:sz(1)
                tpsd = nan(1,sz(2));
                PSDt = squeeze(PSD(it,:,:,:));
                for il=1:sz(2)
                    psdvec = squeeze(PSDt(il,:,:));
                    muvec= squeeze(mu(il,:,:));
                    kvec = squeeze(K(il,:,:));
                    test = griddata(muvec,kvec,psdvec,log10(target_mu),target_K.^(1/3),'natural');
                    tpsd(il) = test;
                end
                PSD_new(it,:) = tpsd;
                if mod(it,100) == 0
                    t2 = toc(t1);
                    ltime = t2/it;
                    est = (sz(1)-it) * ltime;
                    fprintf('%i/%i Estimated Time: %s min\n',it,sz(1),num2str(est/60));
                end
            end
            PSD_new = 10.^PSD_new;
            fprintf('\n');
            fprintf('done\n');
        end    

        psd_rean = PSD_new ./ 2.997e7; %VERB units to (c/cm/MeV).^3

        psd_month = psd_rean;
        if i ~= 1%2
            psd_LT = vertcat(psd_LT,psd_month(1:end-1,:));
        elseif i == 1%2
            psd_LT = vertcat(psd_LT,psd_month(1:end,:));        
        end
    %      Kp_LT = vertcat(Kp_LT,Kp(1:end-1)');

    end


    % Save file
    if m == 1; target_mu_filename = '2000'; end
    if m == 2; target_mu_filename = '3000'; end
    if m == 3; target_mu_filename = '4000'; end
    if m == 4; target_mu_filename = '5000'; end
    if m == 5; target_mu_filename = '6000'; end
    if m == 6; target_mu_filename = '7000'; end
    if m == 7; target_mu_filename = '8000'; end
    if m == 8; target_mu_filename = '9000'; end
    if m == 9; target_mu_filename = '10000'; end
    if m == 10; target_mu_filename = '12000'; end
    if m == 11; target_mu_filename = '15000'; end

%     filename = ['psdfile_rean_',target_mu_filename,'MeVG_0.1GRe_3D+MT+EMICS+MP_Oct12_Chorus2018Hiss2016EMICS2_q=5000r.mat'];
%     filename = ['psdfile_rean_',target_mu_filename,'MeVG_0.01GRe_3D+MT+EMICS+MP_Oct12Sep13.mat'];
    filename = ['psdfile_obs_',target_mu_filename,'MeVG_0.1GRe_3D+MT+EMICS+MP_Oct12_FAST_q=5000r.mat'];
    save(filename,'psd_LT')


end


