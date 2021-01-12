function Res = BCCT_Perm_CaSCNMat_Cal(RealCompPara,ROIsignals)
% ROIsignals: 1*N
Calmethod1 = RealCompPara.Calmethod1;
COVcond = RealCompPara.COVcond;
if COVcond(1)==1
    theCovariables = [ones(size(ROIsignals,1),1),RealCompPara.COVs];
else
    theCovariables = ones(size(ROIsignals,1),1);
end
IOrds = RealCompPara.IMGORD;
Order = RealCompPara.GCAorder;
if Calmethod1 % res
    for i = 1:size(ROIsignals,2)
        [ResultMap1,ResultMap2,ResultMap3,ResultMap4,ResultMap5] = restgca_residual(ROIsignals(IOrds,i),ROIsignals(IOrds,:),Order,theCovariables(IOrds,:));
        ResultMap_p1 = wgr_pwGC_F(ResultMap1,size(ROIsignals,1),Order);
        ResultMap_p2 = wgr_pwGC_F(ResultMap2,size(ROIsignals,1),Order);
        Res.res.ResultMap1(i,:) = ResultMap1;
        Res.res.ResultMap2(i,:) = ResultMap2;
        Res.res.ResultMap3(i,:) = ResultMap3;
        Res.res.ResultMap4(i,:) = ResultMap4;
        Res.res.ResultMap5(i,:) = ResultMap5;
        Res.res.ResultMap_p1(i,:) = ResultMap_p1;
        Res.res.ResultMap_p2(i,:) = ResultMap_p2;
    end
else % coef
    for i = 1:size(ROIsignals,2)
        [ResultMap1,ResultMap2,ResultMap3,ResultMap4,res]=restgca_coefficient_XQ(ROIsignals(IOrds,i), ROIsignals(IOrds,:),Order,theCovariables(IOrds,:));
        for j = 1:Order
            ResMap1(j,i,:) = ResultMap1{j};
            ResMap2(j,i,:) = ResultMap2{j};
            ResT1(j,i,:) = res.T_T1{j};
            ResT2(j,i,:) = res.T_T2{j};
            ResZ1(j,i,:) = res.Z_T1{j};
            ResZ2(j,i,:) = res.Z_T2{j};
            ResP1(j,i,:) = res.P_T1{j};
            ResP2(j,i,:) = res.P_T2{j};
        end
    end
    Res.coef.ResultMap1 = ResMap1;
    Res.coef.ResultMap2 = ResMap2;
    Res.coef.ResT1 = ResT1;
    Res.coef.ResT2 = ResT2;
    Res.coef.ResZ1 = ResZ1;
    Res.coef.ResZ2 = ResZ2;
    Res.coef.ResP1 = ResP1;
    Res.coef.ResP2 = ResP2;
    end
end

if RealCompPara.PermMark
    PermNum = RealCompPara.PermNum;
    clear ResultMap1 ResultMap2 ResultMap3 ResultMap4 ResultMap5 
    for iperm = 1:PermNum
        IOrds = randperm(size(ROIsignals,1));
        if Calmethod1 % res
            for i = 1:size(ROIsignals,2)
                [ResultMap1(iperm,i,:),ResultMap2(iperm,i,:),ResultMap3(iperm,i,:),ResultMap4(iperm,i,:),ResultMap5(iperm,i,:)] = ...
                    restgca_residual(ROIsignals(IOrds,i), ROIsignals(IOrds,:),Order,theCovariables(IOrds,:));
            end
        else
            for i = 1:size(ROIsignals,2)
                [ResultMap1,ResultMap2,ResultMap3,ResultMap4,~]=restgca_coefficient_XQ(ROIsignals(IOrds,i), ROIsignals(IOrds,:),Order,theCovariables(IOrds,:));
                for j = 1:Order
                    ResMap1(j,iperm,i,:) = ResultMap1{j};
                    ResMap2(j,iperm,i,:) = ResultMap2{j};
                end
            end
        end
    end
    if Calmethod1
        for i = 1:size(ROIsignals,2)
            [mu,sig,~,sigci] = normfit(squeeze(ResultMap1(:,i,:)));
            P_map1(i,:) = normcdf(Res.res.ResultMap1(i,:),mu,sig);
            
            [mu,sig,~,sigci] = normfit(squeeze(ResultMap2(:,i,:)));
            P_map2(i,:) = normcdf(Res.res.ResultMap2(i,:),mu,sig);;
            
            [mu,sig,~,sigci] = normfit(squeeze(ResultMap3(:,i,:)));
            P_map3(i,:) = normcdf(Res.res.ResultMap3(i,:),mu,sig);
            
            [mu,sig,~,sigci] = normfit(squeeze(ResultMap4(:,i,:)));
            P_map4(i,:) = normcdf(Res.res.ResultMap4(i,:),mu,sig);
            
            [mu,sig,~,sigci] = normfit(squeeze(ResultMap5(:,i,:)));
            P_map5(i,:) = normcdf(Res.res.ResultMap5(i,:),mu,sig);
        end
    else
        for iord = 1:Order
            ResultMap1t = squeeze(ResMap1(iord,:,:,:));
            ResultMap2t = squeeze(ResMap2(iord,:,:,:));
            for i = 1:size(ROIsignals,2)                
                [mu,sig,~,sigci] = normfit(squeeze(ResultMap1t(:,i,:)));
                P_map1(iord,i,:) = normcdf(squeeze(Res.coef.ResultMap1(iord,i,:))',mu,sig);
                [mu,sig,~,sigci] = normfit(squeeze(ResultMap2t(:,i,:)));
                P_map2(iord,i,:) = normcdf(squeeze(Res.coef.ResultMap2(iord,i,:))',mu,sig);
            end
        end
    end
    if Calmethod1
        Res.res.P_map1 = P_map1;
        Res.res.P_map2 = P_map2;
        Res.res.P_map3 = P_map3;
        Res.res.P_map4 = P_map4;
        Res.res.P_map5 = P_map5;
    else
        Res.coef.P_map1 = P_map1;
        Res.coef.P_map2 = P_map2;
    end
end
end
