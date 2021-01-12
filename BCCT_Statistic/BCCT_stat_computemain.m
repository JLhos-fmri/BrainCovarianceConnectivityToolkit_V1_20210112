function BCCT_stat_computemain(Parameter)
pathfiles = which('BCCT_stat_computemain.m');
[path nam ext] = fileparts(pathfiles);
addpath(fullfile(path,'surfstat'));
Outputdir = Parameter.Outputdir;
RealCompPara.Outputdir = Outputdir;
Inputdir1 = Parameter.Inputdir1;
RealCompPara.Inputdir1 = Inputdir1;
Inputdir2 = Parameter.Inputdir2;
RealCompPara.Inputdir2 = Inputdir2;
Mapmod = Parameter.Mapmod;
RealCompPara.Mapmod = Mapmod;
STTYPE = Parameter.sttype;
RealCompPara.STTYPE = STTYPE;

if Mapmod % mapmod
    maskfiles = Parameter.maskdir;
    [vmask,datamask] = Dynamic_read_dir_NIFTI(maskfiles);
    SIG1dir = fullfile(Inputdir1,'maskedSignal.mat');
    SIG2dir = fullfile(Inputdir2,'maskedSignal.mat');
    ROIsig1dir = fullfile(Inputdir1,'ROIsignal.mat');
    ROIsig2dir = fullfile(Inputdir2,'ROIsignal.mat');
    SIG1 = load(SIG1dir);
    SIG2 = load(SIG2dir);
    roi1 = load(ROIsig1dir);
    roi2 = load(ROIsig2dir);
    SIG1_SIG = SIG1.maskedSignal;
    SIG2_SIG = SIG2.maskedSignal;
    ROI1_SIG = roi1.ROIsignals;
    ROI2_SIG = roi2.ROIsignals;
    MASK0 = SIG1.DATMASK;
    MASK0_T = SIG2.DATMASK;
    if any(MASK0-MASK0_T)
        error('two data sets used different mask');
    end
    indexs = find(MASK0);
    NumG1 = size(SIG1_SIG,1);
    NumG2 = size(SIG2_SIG,1);
    Nroi = size(ROI1_SIG,2);
    if STTYPE % interaction
        disp('something wrong for the interaction')
%         for i = 1:NumG1
%             GroupList{i,1} = 'G1';
%         end
%         for i = 1:NumG2
%             GroupList{i+NumG1,1} = 'G2';
%         end
%         Y = [SIG1_SIG;SIG2_SIG];
%         Group = term(GroupList);
%         for i = 1:Nroi
%             ROIsig = [ROI1_SIG(:,i);ROI2_SIG(:,i)];
%             ROITERM = term(ROIsig);
%             XDesign = 1+Group+ROITERM+Group*ROITERM;
%             slm = SurfStatLinMod(Y,XDesign);
%             %             slm = SurfStatT(slm,ROITERM*Group.G1-ROITERM*Group.G2);
%             slm = SurfStatT(slm,[0 0 0 0 1 -1]);
%             %             pval.P=stat_threshold(0,1,0,df,[10 slm.t(1)],[],[],[],slm.k,[],[],0)
%             [pval] = stat_threshold( 0, 1, 0, slm.df, [10,slm.t],[], [], [],slm.k, [], [], 0 );
%             PVALS = pval(2:end);
%             Tval(i,:) = slm.t;
%             Pval(i,:) = PVALS;
%             clear ROIsig ROITERM XDesign slm pval PVALS
%         end
%         TPvaldir = fullfile(Outputdir,'T_Pval.mat');
%         save(TPvaldir,'Tval','Pval');
%         for i = 1:Nroi
%             if i<10
%                 OutfilenametempT = fullfile(Outputdir,['Tmap_ROI00000',num2str(i),'.nii']);
%                 OutfilenametempP = fullfile(Outputdir,['Pmap_ROI00000',num2str(i),'.nii']);
%             elseif i<100
%                 OutfilenametempT = fullfile(Outputdir,['Tmap_ROI0000',num2str(i),'.nii']);
%                 OutfilenametempP = fullfile(Outputdir,['Pmap_ROI0000',num2str(i),'.nii']);
%             elseif i<1000
%                 OutfilenametempT = fullfile(Outputdir,['Tmap_ROI000',num2str(i),'.nii']);
%                 OutfilenametempP = fullfile(Outputdir,['Pmap_ROI000',num2str(i),'.nii']);
%             else
%                 OutfilenametempT = fullfile(Outputdir,['Tmap_ROI00',num2str(i),'.nii']);
%                 OutfilenametempP = fullfile(Outputdir,['Pmap_ROI00',num2str(i),'.nii']);
%             end
%             DATOUT = zeros(size(MASK0));
%             DATOUT(indexs) = Tval(i,:);
%             DynamicBC_write_NIFTI(DATOUT,vmask,OutfilenametempT);
%             DATOUT = zeros(size(MASK0));
%             DATOUT(indexs) = Pval(i,:);
%             DynamicBC_write_NIFTI(DATOUT,vmask,OutfilenametempP);
%         end
    else % perm
        COVCOND = Parameter.COVCOND;
        if COVCOND
            COV1 = load(Parameter.Cov1dir);
            COV2 = load(Parameter.Cov2dir);
            COVTOTAL = [COV1;COV2];
            [R1 P1] = partialcorr(SIG1_SIG,ROI1_SIG,COV1);
            [R2 P2] = partialcorr(SIG2_SIG,ROI2_SIG,COV2);
            DF_E1 = size(SIG1_SIG,1)-2-size(COV1,2);
            DF_E2 = size(SIG1_SIG,1)-2-size(COV1,2);
        else
            [R1 P1] = corr(SIG1_SIG,ROI1_SIG);
            [R2 P2] = corr(SIG2_SIG,ROI2_SIG);
            DF_E1 = size(SIG1_SIG,1)-2;
            DF_E2 = size(SIG1_SIG,1)-2;
        end
        [Z1, ZP1] = AS_TFRtoZ(R1,'R',DF_E1,[]);
        [Z2, ZP2] = AS_TFRtoZ(R2,'R',DF_E2,[]);
        RDIFF = R1-R2;
        ZDIFF = Z1-Z2;
        TotalNum = NumG1+NumG2;
        permtimes = Parameter.permtime;
        Gsig = [SIG1_SIG;SIG2_SIG];
        Rsig = [ROI1_SIG;ROI2_SIG];
        for i = 1:Nroi
            rdiff = zeros(size(SIG1_SIG,2),permtimes);
            zdiff = zeros(size(SIG1_SIG,2),permtimes);
            for iperm = 1:permtimes
                randord = randperm(TotalNum);
                G1Randp = Gsig(randord(1:NumG1),:);
                G2Randp = Gsig(randord(NumG1+1:TotalNum),:);
                R1Randp = Rsig(randord(1:NumG1),i);
                R2Randp = Rsig(randord(NumG1+1:TotalNum),i);
                if COVCOND
                    COV1Randp = COVTOTAL(randord(1:NumG1),:);
                    COV2Randp = COVTOTAL(randord(NumG1+1:TotalNum),:);
                    [r1 p1] = partialcorr(G1Randp,R1Randp,COV1Randp);
                    [r2 p2] = partialcorr(G2Randp,R2Randp,COV2Randp);
                else
                    [r1 p1] = corr(G1Randp,R1Randp);
                    [r2 p2] = corr(G2Randp,R2Randp);
                end
                rdiff(:,iperm) = r1-r2;
                
                [z1, zP1] = AS_TFRtoZ(r1,'R',DF_E1,[]);
                [z2, zP2] = AS_TFRtoZ(r2,'R',DF_E2,[]);
                zdiff(:,iperm) = z1-z2;
            end
            RDIFFSEP = RDIFF(:,i);
            [mu,sig,~,sigci] = normfit(rdiff');
            P = normcdf(RDIFFSEP',mu,sig);
            
            ZDIFFSEP = ZDIFF(:,i);
            [mu,sig,~,sigci] = normfit(zdiff');
            ZP = normcdf(ZDIFFSEP',mu,sig);
            
            if i<10
                OutfilenametempP = fullfile(Outputdir,['ZPmap_ROI00000',num2str(i),'.nii']);
                OutfilenametempPZ = fullfile(Outputdir,['ZPmap_ROI00000',num2str(i),'.nii']);
            elseif i<100
                OutfilenametempP = fullfile(Outputdir,['ZPmap_ROI0000',num2str(i),'.nii']);
                OutfilenametempPZ = fullfile(Outputdir,['ZPmap_ROI0000',num2str(i),'.nii']);
            elseif i<1000
                OutfilenametempP = fullfile(Outputdir,['ZPmap_ROI000',num2str(i),'.nii']);
                OutfilenametempPZ = fullfile(Outputdir,['ZPmap_ROI000',num2str(i),'.nii']);
            else
                OutfilenametempP = fullfile(Outputdir,['ZPmap_ROI00',num2str(i),'.nii']);
                OutfilenametempPZ = fullfile(Outputdir,['ZPmap_ROI00',num2str(i),'.nii']);
            end
            DATOUT = zeros(size(MASK0));
            DATOUT(indexs) = P;
            DynamicBC_write_NIFTI(DATOUT,vmask,OutfilenametempP);
            DATOUT = zeros(size(MASK0));
            DATOUT(indexs) = ZP;
            DynamicBC_write_NIFTI(DATOUT,vmask,OutfilenametempPZ);
        end
    end
else % matrixmod
    Parlab = load(fullfile(Inputdir1,'SetUpparameter.mat'));
    Partinfo = Parlab.Parameter.Partinfo;
    RealCompPara.Partinfo = Partinfo;
    
    SIG1dir = fullfile(Inputdir1,'ROIsignal.mat');
    SIG2dir = fullfile(Inputdir2,'ROIsignal.mat');
    roi1 = load(SIG1dir);
    roi2 = load(SIG2dir);
    ROI1_SIG = roi1.ROIsignals;
    ROI2_SIG = roi2.ROIsignals;
    NumG1 = size(ROI1_SIG,1);
    NumG2 = size(ROI2_SIG,1);
    Nroi = size(ROI1_SIG,2);
    scatterlab1 = dir(fullfile(Inputdir1,'Scatter_R_Pres.mat'));
    scatterlab2 = dir(fullfile(Inputdir2,'Scatter_R_Pres.mat'));
    if STTYPE % interaction 
        disp('something wrong for the interaction')
%         for i = 1:NumG1
%             GroupList{i,1} = 'G1';
%         end
%         for i = 1:NumG2
%             GroupList{i+NumG1,1} = 'G2';
%         end
%         Y = [ROI1_SIG;ROI2_SIG];
%         Group = term(GroupList);
%         for i = 1:Nroi
%             ROIsig = [ROI1_SIG(:,i);ROI2_SIG(:,i)];
%             ROITERM = term(ROIsig);
%             XDesign = 1+Group+ROITERM+Group*ROITERM;
%             slm = SurfStatLinMod(Y,XDesign);
%             %             slm = SurfStatT(slm,ROITERM*Group.G1-ROITERM*Group.G2);
%             slm = SurfStatT(slm,[0 0 0 0 1 -1]);
%             %             pval.P=stat_threshold(0,1,0,df,[10 slm.t(1)],[],[],[],slm.k,[],[],0)
%             [pval] = stat_threshold( 0, 1, 0, slm.df, [10,slm.t],[], [], [],slm.k, [], [], 0 );
%             PVALS = pval(2:end);
%             Tval(i,:) = slm.t;
%             Pval(i,:) = PVALS;
%             clear ROIsig ROITERM XDesign slm pval PVALS
%         end
%         Tval(1:size(Tval)+1:end) = 0;
%         
%         TPvaldir = fullfile(Outputdir,'T_Pval.mat');
%         save(TPvaldir,'Tval','Pval');
    else % perm
        COVCOND = Parameter.COVCOND;
        if COVCOND
            COV1 = load(Parameter.Cov1dir);
            COV2 = load(Parameter.Cov2dir);
            COVTOTAL = [COV1;COV2];
            if Partinfo                
                R1 = zeros(size(ROI1_SIG,2));
                P1 = ones(size(ROI1_SIG,2));  
                R2 = zeros(size(ROI1_SIG,2));
                P2 = ones(size(ROI1_SIG,2));
                for i = 1:size(ROI1_SIG,2)-1
                    for j = i+1:size(ROI1_SIG,2)
                        sigind = 1:size(ROI1_SIG,2);
                        sigind([i,j]) = [];
                        COVnew = [COV1,ROI1_SIG(:,sigind)];
                        [r p] = partialcorr(ROI1_SIG(:,i),ROI1_SIG(:,j),COVnew);
                        R1(i,j) = r;
                        P1(i,j) = p;
                        R1(j,i) = r;
                        P1(j,i) = p;
                        %                        
                        COVnew = [COV2,ROI2_SIG(:,sigind)];
                        [r p] = partialcorr(ROI2_SIG(:,i),ROI2_SIG(:,j),COVnew);
                        R2(i,j) = r;
                        P2(i,j) = p;
                        R2(j,i) = r;
                        P2(j,i) = p;
                    end
                end
                DF_E1 = size(ROI1_SIG,1)-2-size(COVnew,2);
                DF_E2 = size(ROI2_SIG,1)-2-size(COVnew,2);
                [Z1,PZ1] = AS_TFRtoZ(R1,'R',DF_E1,[]);
                [Z2,PZ2] = AS_TFRtoZ(R2,'R',DF_E2,[]);
            else
                [R1 P1] = partialcorr(ROI1_SIG,ROI1_SIG,COV1);
                [R2 P2] = partialcorr(ROI2_SIG,ROI2_SIG,COV2);                
                DF_E1 = size(ROI1_SIG,1)-2-size(COV1,2);          
                DF_E2 = size(ROI2_SIG,1)-2-size(COV2,2);
                [Z1, PZ1] = AS_TFRtoZ(R1,'R',DF_E1,[]);
                [Z2, PZ2] = AS_TFRtoZ(R2,'R',DF_E2,[]);
            end
        else
            if Partinfo                      
                R1 = zeros(size(ROI1_SIG,2));
                P1 = ones(size(ROI1_SIG,2));
                R2 = zeros(size(ROI2_SIG,2));
                P2 = ones(size(ROI2_SIG,2));
                for i = 1:size(ROI1_SIG,2)-1
                    for j = i+1:size(ROI1_SIG,2)
                        sigind = 1:size(ROI1_SIG,2);
                        sigind([i,j]) = [];
                        COVnew = ROI1_SIG(:,sigind);
                        [r p] = partialcorr(ROI1_SIG(:,i),ROI1_SIG(:,j),COVnew);
                        R1(i,j) = r;
                        P1(i,j) = p;
                        R1(j,i) = r;
                        P1(j,i) = p;
                        %                        
                        COVnew = ROI2_SIG(:,sigind);
                        [r p] = partialcorr(ROI2_SIG(:,i),ROI2_SIG(:,j),COVnew);
                        R2(i,j) = r;
                        P2(i,j) = p;
                        R2(j,i) = r;
                        P2(j,i) = p;                        
                    end
                end
                DF_E1 = size(ROI1_SIG,1)-2-size(COVnew,2);
                DF_E2 = size(ROI2_SIG,1)-2-size(COVnew,2);
                [Z1,PZ1] = AS_TFRtoZ(R1,'R',DF_E1,[]);
                [Z2,PZ2] = AS_TFRtoZ(R2,'R',DF_E2,[]);
            else
                [R1 P1] = corr(ROI1_SIG,ROI1_SIG);
                [R2 P2] = corr(ROI2_SIG,ROI2_SIG);
                DF_E1 = size(ROI1_SIG,1)-2;
                DF_E2 = size(ROI2_SIG,1)-2;
                [Z1,PZ1] = AS_TFRtoZ(R1,'R',DF_E1,[]);
                [Z2,PZ2] = AS_TFRtoZ(R2,'R',DF_E2,[]);
            end
        end
        RDIFF = R1-R2;
        ZDIFF = Z1-Z2;
        TotalNum = NumG1+NumG2;
        permtimes = Parameter.permtime;
        Gsig = [ROI1_SIG;ROI2_SIG];
        
%%      
        try
            rdiff = zeros(size(ROI1_SIG,2),size(ROI1_SIG,2),permtimes);
            for iperm = 1:permtimes
                randord = randperm(TotalNum);
                G1Randp = Gsig(randord(1:NumG1),:);
                G2Randp = Gsig(randord(NumG1+1:TotalNum),:);
                if COVCOND
                    COV1Randp = COVTOTAL(randord(1:NumG1),:);
                    COV2Randp = COVTOTAL(randord(NumG1+1:TotalNum),:);
                    if Partinfo                        
                        for i = 1:size(G1Randp,2)-1
                            for j = 1+1:size(G1Randp,2)
                                covindperm = 1:size(G1Randp,2);
                                covindperm([i,j]) = [];
                                COVnewperm1 = [G1Randp(:,covindperm),COV1Randp];
                                COVnewperm2 = [G2Randp(:,covindperm),COV2Randp];
                                [r1perm p1perm] = partialcorr(G1Randp(:,i),G1Randp(:,j),COVnewperm1);
                                [r2perm p2perm] = partialcorr(G2Randp(:,i),G2Randp(:,j),COVnewperm2);
                                r1(i,j) = r1perm;
                                p1(i,j) = p1perm;
                                r1(j,i) = r1perm;
                                p1(j,i) = p1perm;
                                r2(i,j) = r2perm;
                                p2(i,j) = p2perm;
                                r2(j,i) = r2perm;
                                p2(j,i) = p2perm;
                            end
                        end
                    else
                        [r1 p1] = partialcorr(G1Randp,G1Randp,COV1Randp);
                        [r2 p2] = partialcorr(G2Randp,G2Randp,COV2Randp);
                    end
                else
                    if Partinfo
                        for i = 1:size(G1randp,2)-1
                            for j = 1+1:size(G1randp,2)
                                covindperm = 1:size(G1randp,2);
                                covindperm([i,j]) = [];
                                COVnewperm1 = [G1randp(:,covindperm)];
                                COVnewperm2 = [G2randp(:,covindperm)];
                                [r1perm p1perm] = partialcorr(G1Randp(:,i),G1Randp(:,j),COVnewperm1);
                                [r2perm p2perm] = partialcorr(G2Randp(:,i),G2Randp(:,j),COVnewperm2);
                                r1(i,j) = r1perm;
                                p1(i,j) = p1perm;
                                r1(j,i) = r1perm;
                                p1(j,i) = p1perm;
                                r2(i,j) = r2perm;
                                p2(i,j) = p2perm;
                                r2(j,i) = r2perm;
                                p2(j,i) = p2perm;
                            end
                        end
                    else
                        [r1 p1] = corr(G1Randp,G1Randp);
                        [r2 p2] = corr(G2Randp,G2Randp);
                    end
                end
                [Z1 PZ1] = AS_TFRtoZ(r1,'R',DF_E1,[]); 
                [Z2 PZ2] = AS_TFRtoZ(r2,'R',DF_E2,[]); 
                rdiff(:,:,iperm) = r1-r2;
                zdiff(:,:,iperm) = Z1-Z2;
            end
            for i = 1:Nroi
                RDIFFSEP = RDIFF(:,i);
                rdiff2 = squeeze(rdiff(:,i,:));
                [mu,sig,muci,sigci] = normfit(rdiff2');
                P = normcdf(RDIFFSEP',mu,sig);
                Pout(i,:) = P;               
                
                ZDIFFSEP = ZDIFF(:,i);
                zdiff2 = squeeze(zdiff(:,i,:));
                [mu,sig,muci,sigci] = normfit(zdiff2');
                ZP = normcdf(ZDIFFSEP',mu,sig);
                ZPout(i,:) = ZP;
            end
        catch
            for i = 1:Nroi
                rdiff = zeros(size(ROI1_SIG,2),permtimes);
                zdiff = zeros(size(ROI1_SIG,2),permtimes);
                for iperm = 1:permtimes
                    randord = randperm(TotalNum);
                    G1Randp = Gsig(randord(1:NumG1),:);
                    G2Randp = Gsig(randord(NumG1+1:TotalNum),:);
                    R1Randp = Gsig(randord(1:NumG1),i);
                    R2Randp = Gsig(randord(NumG1+1:TotalNum),i);
                    if COVCOND
                        COV1Randp = COVTOTAL(randord(1:NumG1),:);
                        COV2Randp = COVTOTAL(randord(NumG1+1:TotalNum),:);
                        if Partinfo
                            for j = 1:size(ROI1_SIG,2)
                                indsiguse = 1:size(ROI1_SIG,2);
                                if i~=j
                                    indsiguse([i,j]) = [];
                                else
                                    indsiguse(i) = [];
                                end                                
                                COVperm1orig = [G1Randp(:,indsiguse)];
                                COVperm2orig = [G2Randp(:,indsiguse)];                                
                                COVnewperm1 = [COVperm1orig,COV1Randp];
                                COVnewperm2 = [COVperm2orig,COV2Randp];
                                [r1o p1o] = partialcorr(G1Randp(:,j),R1Randp,COVnewperm1);
                                [r2o p2o] = partialcorr(G2Randp(:,j),R2Randp,COVnewperm2);
                                r1(j) = r1o;
                                r2(j) = r2o;
                            end                                
                        else
                            [r1 p1] = partialcorr(G1Randp,R1Randp,COV1Randp);
                            [r2 p2] = partialcorr(G2Randp,R2Randp,COV2Randp);
                        end
                    else
                        if Partinfo
                            for j = 1:size(ROI1_SIG,2)
                                indsiguse = 1:size(ROI1_SIG,2);
                                if i~=j
                                    indsiguse([i,j]) = [];
                                else
                                    indsiguse(i) = [];
                                end                                
                                COVperm1orig = [G1Randp(:,indsiguse)];
                                COVperm2orig = [G2Randp(:,indsiguse)];                                
                                COVnewperm1 = [COVperm1orig];
                                COVnewperm2 = [COVperm2orig];
                                [r1o p1o] = partialcorr(G1Randp(:,j),R1Randp,COVnewperm1);
                                [r2o p2o] = partialcorr(G2Randp(:,j),R2Randp,COVnewperm2);
                                r1(j) = r1o;
                                r2(j) = r2o;
                            end
                        else
                            [r1 p1] = corr(G1Randp,R1Randp);
                            [r2 p2] = corr(G2Randp,R2Randp);
                        end
                    end
                    r1(isnan(r1)) = 0;
                    r2(isnan(r2)) = 0;
                    
                    [Z1 PZ1] = AS_TFRtoZ(r1,'R',DF_E1,[]);
                    [Z2 PZ2] = AS_TFRtoZ(r2,'R',DF_E2,[]);
                    rdiff(:,iperm) = r1-r2;
                    zdiff(:,iperm) = Z1-Z2;
                end
                RDIFFSEP = RDIFF(:,i);                
                [mu,sig,muci,sigci] = normfit(rdiff');
                P = normcdf(RDIFFSEP',mu,sig);
                Pout(i,:) = P;
                
                
                ZDIFFSEP = ZDIFF(:,i);                
                [mu,sig,muci,sigci] = normfit(zdiff');
                ZP = normcdf(ZDIFFSEP',mu,sig);
                ZPout(i,:) = ZP;
            end
        end
%%
        Pval = Pout;
        ZPval = ZPout;
        TPvaldir = fullfile(Outputdir,'PermPval.mat');
        save(TPvaldir,'Pval','ZPval');
    end
%     if ~isempty(scatterlab1)&&~isempty(scatterlab2)
%         SCAT1 = load(fullfile(Inputdir1,'Scatter_R_Pres.mat'));
%         SCAT2 = load(fullfile(Inputdir2,'Scatter_R_Pres.mat'));
%         if STTYPE % interaction
%             Tval = zeros(SCAT1.nrois);
%             Pval = zeros(SCAT1.nrois);
%             kind = 1;
%             for i = 1:SCAT1.nrois-1
%                 for j = i+1:SCAT1.nrois
%                     SelLab1 = SCAT1.Sellab(:,kind);
%                     SelLab2 = SCAT2.Sellab(:,kind);
%                     Rsigtemp1_1 = ROI1_SIG(:,i);
%                     Rsigtemp1_2 = ROI1_SIG(:,j);
%                     Rsigtemp2_1 = ROI2_SIG(:,i);
%                     Rsigtemp2_2 = ROI2_SIG(:,j);
%                     RsigX = [Rsigtemp1_1(~SelLab1);Rsigtemp2_1(~SelLab2)];
%                     RsigY = [Rsigtemp1_2(~SelLab1);Rsigtemp2_2(~SelLab2)];
%                     for ix = 1:nnz(~SelLab1)
%                         GroupList{ix,1} = 'G1';
%                     end
%                     for ix = 1:nnz(~SelLab2)
%                         GroupList{ix+nnz(~SelLab1),1} = 'G2';
%                     end
%                     Group = term(GroupList);
%                     Y = RsigY;
%                     ROIsig = RsigX;
%                     ROITERM = term(ROIsig);
%                     XDesign = 1+Group+ROITERM+Group*ROITERM;
%                     slm = SurfStatLinMod(Y,XDesign);
%                     slm = SurfStatT(slm,[0 0 0 0 1 -1]);
%                     [pval] = stat_threshold( 0, 1, 0, slm.df, [10,slm.t],[], [], [],slm.k, [], [], 0 );
%                     PVALS = pval(2:end);
%                     Tval(i,j) = slm.t;
%                     Pval(i,j) = PVALS;
%                     clear ROIsig ROITERM XDesign slm pval PVALS GroupList Group
%                     kind = kind+1;
%                 end
%             end
%             TPvaldir = fullfile(Outputdir,'Scatter_T_Pval.mat');
%             save(TPvaldir,'Tval','Pval');
%         else % perm
%             
%             %             COVCOND = Parameter.COVCOND;
%             COVCOND = 0;% current version, no cov;
%             if COVCOND
%                 COV1 = load(Parameter.Cov1dir);
%                 COV2 = load(Parameter.Cov2dir);
%                 COVTOTAL = [COV1;COV2];
%                 [R1 P1] = partialcorr(ROI1_SIG,ROI1_SIG,COV1);
%                 [R2 P2] = partialcorr(ROI2_SIG,ROI2_SIG,COV2);
%             else
%                 R1 = SCAT1.R;
%                 R2 = SCAT2.R;
%             end
%             RDIFF = R1-R2;
%             permtimes = Parameter.permtime;
%             
%             Pout = ones(SCAT1.nrois)*0.5;
%             kind = 1;
%             for i = 1:SCAT1.nrois-1
%                 for j = i+1:SCAT1.nrois
%                     clear rdiff
%                     for iperm = 1:permtimes
%                         SelLab1 = SCAT1.Sellab(:,kind);
%                         SelLab2 = SCAT2.Sellab(:,kind);
%                         ROISIGT1_1 = ROI1_SIG(~SelLab1,i);
%                         ROISIGT1_2 = ROI1_SIG(~SelLab1,j);
%                         ROISIGT2_1 = ROI2_SIG(~SelLab2,i);
%                         ROISIGT2_2 = ROI2_SIG(~SelLab2,j);
%                         Nsig1 = length(ROISIGT1_1);
%                         Nsig2 = length(ROISIGT2_1);
%                         TotalNum = Nsig1+Nsig2;
%                         randorder = randperm(TotalNum);
%                         ROISIGT1 = [ROISIGT1_1;ROISIGT2_1];
%                         ROISIGT2 = [ROISIGT1_2;ROISIGT2_2];
%                         
%                         [r1 p1] = corr(ROISIGT1(randorder(1:Nsig1)),ROISIGT2(randorder(1:Nsig1)));
%                         [r2 p2] = corr(ROISIGT1(randorder(Nsig1+1:TotalNum)),ROISIGT2(randorder(Nsig1+1:TotalNum)));
%                         rdiff(:,iperm) = r1-r2;
%                     end
%                     [mu,sig,muci,sigci] = normfit(rdiff');
%                     P = normcdf(RDIFF(i,j),mu,sig);
%                     Pout(i,j) = P;
%                     Pout(j,i) = P;
%                     clear rdiff mu sig muci sigci rdiff;
%                     kind = kind+1;
%                 end
%             end
%             
%             Pval = Pout;
%             TPvaldir = fullfile(Outputdir,'Scatter_PermPval.mat');
%             save(TPvaldir,'Pval');
%         end
%     end
end
disp('Compute Finished!')
end