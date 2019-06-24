function [dataOut_C, dataOut_CT, dataOut_CT2]= prepareDataOut(dataIn,cellBody_L_Complete,cellNuclei,cellTubules_L)

%% Just the data and cell/nuclei
dataOut_C                       = dataIn;
% delineate cell and nuclei
dilatedCellNuc                  = uint8(imdilate( zerocross(cellNuclei-(cellBody_L_Complete>0)),ones(2)));
dataOut_C                       = dataOut_C.*(repmat(1-dilatedCellNuc,[1 1 3]))+repmat(255*dilatedCellNuc,[1 1 3]);
%imagesc(dataOut_C)


%% add tubules all with the same colour
numTubules  = max(cellTubules_L(:));
dataOut_CT                     = dataOut_C;
for counterTub                  = 1:numTubules
    dataOut_CT(cellTubules_L==counterTub) =255;
end

%% add the tubules with colours per cells
%clear jet3
jet3                            = zeros(numTubules,3);
jet3(1:end,:)                   = 0.2+round(-15+120*rand(numTubules,3))/100;
jet3(jet3>1)                    = 1;

%%
dataOut_CT2                      = dataOut_C;
for counterTub =1:numTubules
    dataOut_CT2(:,:,1)           = dataOut_CT2(:,:,1).*uint8(cellTubules_L~=counterTub) + uint8(cellTubules_L==counterTub)*jet3(counterTub,1)*255;
    dataOut_CT2(:,:,2)           = dataOut_CT2(:,:,2).*uint8(cellTubules_L~=counterTub) + uint8(cellTubules_L==counterTub)*jet3(counterTub,2)*255;
    dataOut_CT2(:,:,3)           = dataOut_CT2(:,:,3).*uint8(cellTubules_L~=counterTub) + uint8(cellTubules_L==counterTub)*jet3(counterTub,3)*255;
end
%imagesc(dataOut_CT)

%% Tubules with colours per class: cell/no cell/connect
dataOut_CT3                      = dataOut_C;

