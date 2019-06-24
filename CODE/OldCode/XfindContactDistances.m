imagesc(a(:,:,:,1))

%%

imagesc((1-finalCells(:,:,3)).*double(dataIn(:,:,2)))

%%
imagesc (distFromCell(:,:,3))

%%
k1=2;
k2=1;
imagesc((assignedNuclei(:,:,1)==k1).*(distFromCell(:,:,k2))+50*finalCells(:,:,k2)); colorbar
%imagesc((finalCells(:,:,1)).*(distFromCell(:,:,2))); colorbar

%%
D_CT = zeros(numNuclei);
D_CC = zeros(numNuclei);
%%
thresholdDistance                   =50;
for k1=1:numNuclei;
    for k2=1:numNuclei;
        if k1~=k2
            %
            k1=1;k2=2;
            currentDistances                              = distFromCell(:,:,k1);
            distOverTubules                 = currentDistances.*(assignedNuclei==k2);
            connectingTubulesIm             = microTubules_L((distOverTubules<thresholdDistance)&(distOverTubules>0));
            connectingTubules               = unique(connectingTubulesIm);
            [yd1,xd1]=hist(currentDistances(assignedNuclei==k2),(1:200));
            [yd2,xd2]=hist(currentDistances(finalCells(:,:,k2)>0),(1:200));
            %[yd3,xd3]=hist(currentDistances(assignedNuclei==2),(1:300));
            plot(xd1,(yd1)/sum(yd1),'b-',xd2,(yd2)/sum(yd2),'r--')
            grid on
            D_CT(k1,k2) = sum(yd1(1:thresholdDistance))/sum(yd1);
            D_CC(k1,k2) = sum(yd2(1:thresholdDistance))/sum(yd2);
            %
        end
    end
end
