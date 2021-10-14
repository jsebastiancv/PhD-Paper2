function [new_arr] = fast_interpolate_noL(epc_arr, alpha_arr, f_arr, new_epc, new_alpha, varargin)
    loginterp = false;
    if (size(varargin) > 0) 
        if (strcmp(varargin(1), 'log'))
            loginterp = true;
            f_arr = log10(f_arr);
%            for it = 1:size(f_arr, 1)
%                f_arr(it, :, :, :) = log10(f_arr(it, :, :, :));
%            end
        end
    end

    for it = 1:size(f_arr, 1)
        
% 		if (rem(it,100) == 0)
            %system('tput cuu 1');
% 			fprintf('%s/%s\n',['interpolation step = ', num2str(it)],...
%            num2str(size(f_arr,1)))            
% 		end
% 	
        if (loginterp) 
           new_arr(it, 1) = -22; %Trick - PSD=0 for L=1
        else
           new_arr(it, 1) = 1e-22; %Trick - PSD=0 for L=1
        end
        for L_ind = 2:size(f_arr, 2)

            try 
                
            % using the fact, that alpha is independent of epc!
            alpha_l_ind = max(find(alpha_arr(L_ind,1,:)<new_alpha(L_ind)));
            alpha_r_ind = min(find(alpha_arr(L_ind,1,:)>new_alpha(L_ind)));

            % finding bottom and top indexes of epc for each alpha
            epc_lb_ind = max(find(epc_arr(L_ind,:,alpha_l_ind)<new_epc(L_ind)));
            epc_lt_ind = min(find(epc_arr(L_ind,:,alpha_l_ind)>new_epc(L_ind)));
            epc_rb_ind = max(find(epc_arr(L_ind,:,alpha_r_ind)<new_epc(L_ind)));
            epc_rt_ind = min(find(epc_arr(L_ind,:,alpha_r_ind)>new_epc(L_ind)));

            %['L=', num2str(L_arr(L_ind,1,1)), '; al=', num2str(alpha_l_ind), '; ar=', num2str(alpha_r_ind),'; epc_lb=', num2str(epc_lb_ind), '; epc_lt=', num2str(epc_lt_ind), '; epc_rb=', num2str(epc_rb_ind), '; epc_rt=', num2str(epc_rt_ind)]
            
            % interpolating
%				double a = (x - m_x[i - 1]) / (m_x[i] - m_x[i - 1]);
%				return m_y[i - 1] + a * (m_y[i] - m_y[i - 1]);
            % left side
            fb = f_arr(it, L_ind, epc_lb_ind, alpha_l_ind);
            ft = f_arr(it, L_ind, epc_lt_ind, alpha_l_ind);
            epcb = epc_arr(L_ind, epc_lb_ind, alpha_l_ind);
            epct = epc_arr(L_ind, epc_lt_ind, alpha_l_ind);
            
            a = ( new_epc(L_ind) - epcb ) / ( epct - epcb );
            new_arr_l = fb + a * ( ft - fb );

            % right side
            fb = f_arr(it, L_ind, epc_rb_ind, alpha_r_ind);
            ft = f_arr(it, L_ind, epc_rt_ind, alpha_r_ind);
            epcb = epc_arr(L_ind, epc_rb_ind, alpha_r_ind);
            epct = epc_arr(L_ind, epc_rt_ind, alpha_r_ind);
            
            a = ( new_epc(L_ind) - epcb ) / ( epct - epcb );
            new_arr_r = fb + a * ( ft - fb );

            % between sides
            alphab = alpha_arr(L_ind, epc_lb_ind, alpha_l_ind);
            alphat = alpha_arr(L_ind, epc_rb_ind, alpha_r_ind);
            
            a = ( new_alpha(L_ind) - alphab ) / ( alphat - alphab );
            new_arr(it, L_ind) = new_arr_l + a * ( new_arr_r - new_arr_l );
            
            % temp!!!
            %new_arr(it, L_ind) = f_arr(it, L_ind, epc_rt_ind, alpha_r_ind);
            catch
                err = lasterror;
   %             [err.message, '; L_ind=', num2str(L_ind)]
                if (loginterp) 
                   new_arr(it, L_ind) = -22; %Trick - PSD=0 for L=1
                else
                   new_arr(it, L_ind) = 1e-22; %Trick - PSD=0 for L=1
                end
            end
        end            
        
    end

    if (loginterp) 
        new_arr = 10.^new_arr;
    end
    
end
