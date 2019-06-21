function [dataOut_C, dataOut_CT]= prepareDataOut(dataIn,cellBody_L_Complete,cellNuclei,cellTubules_L)




dataOut_C                   = dataIn;
% delineate cell and nuclei
dilatedCellNuc              = uint8(imdilate( zerocross(cellNuclei-(cellBody_L_Complete>0)),ones(2)));
dataOut_C                   = dataOut_C.*(repmat(1-dilatedCellNuc,[1 1 3]))+repmat(255*dilatedCellNuc,[1 1 3]);
imagesc(dataOut_C)
%% add the tubules with colours
numTubules  = max(cellTubules_L(:));
clear jet3
jet3 = zeros(numTubules,3);
        jet3(1:end,:)=0.2+round(-15+120*rand(numTubules,3))/100;
        jet3(jet3>1)=1;
%%
dataOut_CT                  = dataOut_C;
for counterTub =1:numTubules
    dataOut_CT(:,:,1)           = dataOut_CT(:,:,1).*uint8(cellTubules_L~=counterTub) + uint8(cellTubules_L==counterTub)*jet3(counterTub,1)*255;
    dataOut_CT(:,:,2)           = dataOut_CT(:,:,2).*uint8(cellTubules_L~=counterTub) + uint8(cellTubules_L==counterTub)*jet3(counterTub,2)*255;
    dataOut_CT(:,:,3)           = dataOut_CT(:,:,3).*uint8(cellTubules_L~=counterTub) + uint8(cellTubules_L==counterTub)*jet3(counterTub,3)*255;
end
imagesc(dataOut_CT)


