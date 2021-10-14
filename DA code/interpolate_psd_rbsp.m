
function ypsd = interpolate_psd_rbsp(psd, psd_mu, psd_k, targ_mu, targ_k)

         %we have the rbsp psd at a number of values of mu and k which
         %change with time

         %make array to hold the new psd
         %to start, we will make it such that only one value of mu and k
         %can be requested. - change this later

         ypsd = zeros(length(psd(:, 1, 1)), 1);

         ntime = length(psd(:, 1, 1));
         nmu   = length(psd_mu(1, :, 1));
         nk    = length(psd_k(1, :));

         for tt = 1:ntime

                 %at each time pull out the k array
                 tempk = squeeze(psd_k(tt, :));

                 %check whether the target k is gtr than or less than the
                 %max and min values - if it is, then PSD must be nan
                 if (min(tempk) > targ_k) || (max(tempk) < targ_k)
                    ypsd(tt) = NaN;
                 else
                    tempmu1 = squeeze(psd_mu(tt, :, :));
                    tempmu2 = zeros(nmu, 1);

                    temp_psd1 = squeeze(psd(tt, :, :));
                    temp_psd2 = zeros(nmu, 1);

                    N = ~isnan(tempk);

                    tempk = tempk(N);
                    tempmu1 = tempmu1(:, N);
                    temp_psd1 = temp_psd1(:, N);

                    for tnn = 1:nmu
                        tempmu2(tnn)   = interp1(tempk, squeeze(tempmu1(tnn, :)), targ_k);
                        temp_psd2(tnn) = interp1(tempk, squeeze(temp_psd1(tnn, :)), targ_k);
                    end

                     N = ~isnan(tempmu2);

                     tempmu2 = tempmu2(N);
                     temp_psd2 = temp_psd2(N);

                     if ~isempty(tempmu2)


                            if (min(tempmu2) > targ_mu) || (max(tempmu2) < targ_mu)
                               ypsd(tt) = NaN;
                            else
                               ypsd(tt) = interp1(tempmu2, temp_psd2, targ_mu);
                            end
                     else
                         ypsd(tt) = NaN;
                     end
                 end

         end

end

