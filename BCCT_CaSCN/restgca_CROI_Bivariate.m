function [Result_X2Y,Result_Y2X,ROI_sequence] = restgca_CROI_Bivariate(theROITimeCourses,Order,theCovariables)
nDim4=size(theROITimeCourses,1);
numROIs=size(theROITimeCourses,2);
ROI_sequence=combntns(1:numROIs,2);
Past_1=zeros(nDim4-Order,Order);
Past_2=zeros(nDim4-Order,Order);
Result_X2Y=zeros(Order*2,size(ROI_sequence,1))';
Result_Y2X=zeros(Order*2,size(ROI_sequence,1))';
theCovariables=[theCovariables(Order+1:end,:),ones(nDim4-Order,1)];
for i=1:size(ROI_sequence,1),
    ROI_used=theROITimeCourses(:,ROI_sequence(i,:));
    Now=ROI_used(Order+1:end,:);
    for j=1:Order,
        Past_1(:,j)=ROI_used(j:nDim4-Order+j-1,1);
        Past_2(:,j)=ROI_used(j:nDim4-Order+j-1,2);
        Regressors1=[Past_1,Past_2,theCovariables];
        Regressors2=[Past_2,Past_1,theCovariables];
    end
    b_1=rest_regress(Now(:,2),Regressors1);
    b_2=rest_regress(Now(:,1),Regressors2);
    Result_X2Y(i,:)=b_1(1:Order*2);
    Result_Y2X(i,:)=b_2(1:Order*2);
end
end

function beta = rest_regress(y,X)
         [n,ncolX] = size(X);
         [Q,R,perm] = qr(X,0);
         p = sum(abs(diag(R)) > max(n,ncolX)*eps(R(1)));
         if p < ncolX,
            R = R(1:p,1:p);
            Q = Q(:,1:p);
            perm = perm(1:p);
         end
         beta = zeros(ncolX,1);
         beta(perm) = R \ (Q'*y);
end