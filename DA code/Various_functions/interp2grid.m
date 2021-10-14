function [PSD] = interp2grid(mtime,mInvMu,mInvK,mL,stime,sInvMu,sInvK,sL,sPSD);


nech = size(sInvMu,2);
lvals = mL(:,end,end);
NL = size(mL,1);
NE = size(mL,2);
NA = size(mL,3);
NT = length(mtime);

PSD=nan(NT,NL,NE,NA);
d_mtime = mtime(2) - mtime(1);
dl = lvals(2)-lvals(1);
fprintf(' \n');
parfor it=1:NT
    loop_counter(it,NT,'tput');
    tinx = find(stime >= mtime(it)-d_mtime/2 & stime < mtime(it) + d_mtime/2);
    if ~isempty(tinx)
        spc_InvMu = reshape(sInvMu(tinx,:,:),[],1);
        spc_PSD = reshape(sPSD(tinx,:,:),[],1);

        spc_Lstar = nan(length(tinx),nech,size(sL,2));
        spc_InvK = nan(length(tinx),nech,size(sInvK,2));
        for ie=1:nech
            spc_Lstar(:,ie,:) = sL(tinx,:);
            spc_InvK(:,ie,:) = sInvK(tinx,:);
        end     
        spc_Lstar = reshape(spc_Lstar,[],1);
        spc_InvK = reshape(spc_InvK,[],1);

        nn = find(~isnan(spc_InvMu) & ~isnan(spc_InvK) & ~isnan(spc_Lstar)...
            & ~isnan(spc_PSD));
        if ~isempty(nn)
            spc_Lstart = spc_Lstar(nn);
            spc_InvKt = spc_InvK(nn);
            spc_InvMut = spc_InvMu(nn);
            spc_PSDt = spc_PSD(nn);

            for il=1:NL
                mod_InvMu = reshape(mInvMu(il,:,:),[],1);
                mod_InvK = reshape(mInvK(il,:,:),[],1);
                l1=lvals(il)-dl;
                l2=lvals(il)+dl;
                inx = find(spc_Lstart >= l1 & spc_Lstart < l2);
                if ~isempty(inx)
                    Lt = spc_Lstart(inx);
                    Kt = spc_InvKt(inx);
                    Mut = spc_InvMut(inx);
                    PSDt = spc_PSDt(inx);

                    PSD_out = griddata(log10(Mut),Kt.^(1./3),log10(PSDt),...
                    log10(mod_InvMu), mod_InvK.^(1/3),'natural');
                    if ~isempty(PSD_out)
                        PSD_out = reshape(10.^PSD_out,1,1,...
                            NE,NA);
                        PSD(it,il,:,:) = PSD_out;
                    end
                end
            end 
        end
    end
end
