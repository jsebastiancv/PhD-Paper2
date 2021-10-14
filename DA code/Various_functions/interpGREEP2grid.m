function [PSD] = interp2grid(mtime,mInvMu,mInvK,mL,stime,sInvMu,sInvK,sL,sPSD);

nech = length(sInvMu);
nalph = length(sInvK);
lvals = mL(:,end,end);
NL = size(mL,1);
NE = size(mL,2);
NA = size(mL,3);
NT = length(mtime);

PSD=nan(NT,NL,NE,NA);
d_mtime = mtime(2) - mtime(1);
fprintf(' \n');
parfor it=1:NT
    PSDtmp=nan(NL,NE,NA);
    loop_counter(it,NT,'tput');
    tinx = find(stime >= mtime(it)-d_mtime/2 & stime < mtime(it) + d_mtime/2);
    if ~isempty(tinx)
%        spc_InvMu = reshape(sInvMu(tinx,:,:),[],1);
       spc_PSD = reshape(sPSD(tinx,:,:),[],1);

       spc_InvK = nan(length(tinx),nech,nalph);
       spc_InvMu = nan(length(tinx),nech,nalph);
        for ie=1:nech
            for ia=1:nalph
                spc_InvK(:,ie,ia) = sInvK(ia);
                spc_InvMu(:,ie,ia) = sInvMu(ie);
            end
        end     
        spc_InvK = reshape(spc_InvK,[],1);
        spc_InvMu = reshape(spc_InvMu,[],1);

        nn = find(~isnan(spc_InvMu) & ~isnan(spc_InvK) & ~isnan(spc_PSD));
        if ~isempty(nn)
            spc_InvKt = spc_InvK(nn);
            spc_InvMut = spc_InvMu(nn);
            spc_PSDt = spc_PSD(nn);

            [~,mininx] = min(abs(lvals - sL))

            mod_InvMu = reshape(mInvMu(mininx,:,:),[],1);
            mod_InvK = reshape(mInvK(mininx,:,:),[],1);
            Kt = spc_InvKt;
            Mut = spc_InvMut;
            PSDt = spc_PSDt;

            PSD_out = griddata(log10(Mut),Kt.^(1./3),log10(PSDt),...
            log10(mod_InvMu), mod_InvK.^(1/3),'natural');
            if ~isempty(PSD_out)
                PSD_out = reshape(10.^PSD_out,1,1,...
                    NE,NA);
                PSDtmp(mininx,:,:) = PSD_out;
            end
        end 
    end
    PSD(it,:,:,:) = PSDtmp;
end
