
PSD_BC_temp = squeeze(PSD_Lupper(it_tot,:,:))';
E_VERB_temp = E_VERB';
A_VERB_temp = A_VERB';

t_zone = 0;

col_energy = E_VERB_temp(:);
col_alpha = A_VERB_temp(:);
col_psd = PSD_BC_temp(:);
 
for idx_psd = 1 : length(col_psd)
    
    if isnan(col_psd(idx_psd))
        col_psd(idx_psd) = 0;
    end
end

col_cat = [col_energy,col_alpha,col_psd];

fileID = fopen([simulation_dir,'Input/L_upper_BC.dat'],'w');
    
fprintf(fileID,'VARIABLES = "energy", "alpha", "PSD" \n');

fprintf(fileID,'ZONE T="%f", I=91, J=101 \n',t_zone);

fprintf(fileID,'%e \t %e \t %e \n',col_cat');

fclose(fileID);

