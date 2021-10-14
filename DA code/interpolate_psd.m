
function psd_interp = interpolate_psd(spc_psd, spc_mu, spc_k, target_mu, target_k)
    
    psd_interp = zeros(length(spc_psd(:, 1, 1)), 1);

    ntime = length(spc_psd(:,1,1));
    nmu = length(spc_mu(1,:,1));
%     nk = length(spc_k(1,:));

    for tt = 1 : ntime
        temp_k = squeeze(spc_k(tt,:));  % at each time retrieve the InvK array

        % check whether the target K is greater than or less than the
        % maximum and minimum values, or if InvK is a vector of NaN - if it is, then PSD must be NaN
        if (min(temp_k) > target_k) || (max(temp_k) < target_k || sum(isnan(temp_k)) == length(temp_k))
            psd_interp(tt) = NaN;
        else
            
            temp_mu1 = squeeze(spc_mu(tt,:,:)); % at each time retrieve the InvMuReal array
            temp_mu2 = zeros(nmu, 1);

            temp_psd1 = squeeze(spc_psd(tt,:,:)); % at each time retrieve the PSD array
            temp_psd2 = zeros(nmu, 1);

            N = ~isnan(temp_k);

            temp_k = temp_k(N);
            temp_mu1 = temp_mu1(:,N);
            temp_psd1 = temp_psd1(:,N);

            for tnn = 1:nmu
                temp_mu2(tnn) = interp1(temp_k, squeeze(temp_mu1(tnn,:)), target_k);
                temp_psd2(tnn) = interp1(temp_k, squeeze(temp_psd1(tnn,:)), target_k);
            end

            N = ~isnan(temp_mu2);

            temp_mu2 = temp_mu2(N);
            temp_psd2 = temp_psd2(N);

            if ~isempty(temp_mu2)
                if (min(temp_mu2) > target_mu) || (max(temp_mu2) < target_mu)
                    psd_interp(tt) = NaN;
                else
                    psd_interp(tt) = interp1(temp_mu2, temp_psd2, target_mu);
                end
            else
                psd_interp(tt) = NaN;
            end
        end
    end
end

