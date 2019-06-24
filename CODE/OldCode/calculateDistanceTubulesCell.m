function [D_CT,D_CC,D_TUBULES,dataOut,dataOut2]=calculateDistanceTubulesCell(numNuclei,assignedNuclei,finalCells,distFromCell,microTubules_L,distanceThreshold,dataIn)

%%

D_CT                                        = zeros(numNuclei);
D_CC                                        = zeros(numNuclei);
D_TUBULES{numNuclei,numNuclei}              = 0;

for k1=1:numNuclei;
    %qq                                      = distFromCell(:,:,k1);
    currentDistances                        = distFromCell(:,:,k1);
    for k2=1:numNuclei;
        if k1~=k2
            
            
            distOverTubules                 = currentDistances.*(assignedNuclei==k2);
            connectingTubulesIm             = microTubules_L((distOverTubules<distanceThreshold)&(distOverTubules>0));
            connectingTubules               = unique(connectingTubulesIm);
            [yd1,xd1]                       = hist(currentDistances(assignedNuclei==k2),(1:200));
            [yd2,xd2]                       = hist(currentDistances(finalCells(:,:,k2)>0),(1:200));
            
            %plot(xd1,(yd1)/sum(yd1),'b-',xd2,(yd2)/sum(yd2),'r--')
            %grid on
            D_CT(k1,k2)                     = sum(yd1(1:distanceThreshold))/sum(yd1);
            D_CC(k1,k2)                     = sum(yd2(1:distanceThreshold))/sum(yd2);
            D_TUBULES{k1,k2}                = connectingTubules;
            
        end
    end
end
D_CT = D_CT .*(1-(D_CC>0));


    %% create output with tubules in colours according to nuclei
    
    dataOut                                     = dataIn;
    
        dataOut2                                    = double(dataIn);
    %% Detect direct contact between the cells
    [c1,c2]     = find(D_CC);
    if ~isempty(c1)
        connectingCells                         = sum(finalCells(:,:,c1),3);
        connectingCells                         = (imdilate(zerocross(connectingCells),ones(3)));
        
        dataOut2(:,:,3)                         = dataOut2(:,:,3)+(220*connectingCells);
        %imagesc(uint8(dataOut2))
    end
    %% Detect the tubule that creates the contact between the cells
    [c1,c2]     = find(D_CT);

    if ~isempty(c1)
        for counterContact = 1:numel(c1)
            %firstCell                                   = c1(counterContact);
            %secondCell                                  = c2(counterContact);
            contactTubules                              = D_TUBULES{c1(counterContact),c2(counterContact)}';
            currentTubules                              = ismember(microTubules_L,contactTubules);
            switch (rem(counterContact,7))
                case 0
                    tempdata                            = dataOut2(:,:,1);
                    tempdata(currentTubules)            = 255;
                    dataOut2(:,:,1)                     = tempdata;
                case 1
                    tempdata                            = dataOut2(:,:,2);
                    tempdata(currentTubules)            = 255;
                    dataOut2(:,:,2)                     = tempdata;
                case 2
                    tempdata                            = dataOut2(:,:,3);
                    tempdata(currentTubules)            = 255;
                    dataOut2(:,:,3)                     = tempdata;
                case 3
                    tempdata1                           = dataOut2(:,:,1);
                    tempdata2                           = dataOut2(:,:,2);
                    tempdata1(currentTubules)           = 255;
                    tempdata2(currentTubules)           = 255;
                    dataOut2(:,:,1)                     = tempdata1;
                    dataOut2(:,:,2)                     = tempdata2;
                case 4
                    tempdata1                           = dataOut2(:,:,1);
                    tempdata2                           = dataOut2(:,:,3);
                    tempdata1(currentTubules)           = 255;
                    tempdata2(currentTubules)           = 255;
                    dataOut2(:,:,1)                     = tempdata1;
                    dataOut2(:,:,3)                     = tempdata2;
                    
                case 5
                    tempdata1                           = dataOut2(:,:,3);
                    tempdata2                           = dataOut2(:,:,2);
                    tempdata1(currentTubules)           = 255;
                    tempdata2(currentTubules)           = 255;
                    dataOut2(:,:,3)                     = tempdata1;
                    dataOut2(:,:,2)                     = tempdata2;
                    
                case 6
                    tempdata1                           = dataOut2(:,:,1);
                    tempdata2                           = dataOut2(:,:,2);
                    tempdata3                           = dataOut2(:,:,3);
                    tempdata1(currentTubules)           = 255;
                    tempdata2(currentTubules)           = 255;
                    tempdata3(currentTubules)           = 255;
                    dataOut2(:,:,1)                     = tempdata1;
                    dataOut2(:,:,2)                     = tempdata2;
                    dataOut2(:,:,3)                     = tempdata3;
            end
            
            
            
            
        end
        
        
    end
    
    %imagesc(dataOut2/255)
    
    %% create output with tubules in colours according to nuclei
    
    
    for counterN=1:numNuclei
        switch (rem(counterN,7))
            case 0
                tempdata                            = dataOut(:,:,1);
                tempdata(assignedNuclei==counterN)  = 255;
                dataOut(:,:,1)                      = tempdata;
            case 1
                tempdata                            = dataOut(:,:,2);
                tempdata(assignedNuclei==counterN)  = 255;
                dataOut(:,:,2)                      = tempdata;
            case 2
                tempdata                            = dataOut(:,:,3);
                tempdata(assignedNuclei==counterN)  = 255;
                dataOut(:,:,3)                      = tempdata;
            case 3
                tempdata1                            = dataOut(:,:,1);
                tempdata2                            = dataOut(:,:,2);
                tempdata1(assignedNuclei==counterN)  = 255;
                tempdata2(assignedNuclei==counterN)  = 255;
                dataOut(:,:,1)                      = tempdata1;
                dataOut(:,:,2)                      = tempdata2;
            case 4
                tempdata1                            = dataOut(:,:,1);
                tempdata2                            = dataOut(:,:,3);
                tempdata1(assignedNuclei==counterN)  = 255;
                tempdata2(assignedNuclei==counterN)  = 255;
                dataOut(:,:,1)                      = tempdata1;
                dataOut(:,:,3)                      = tempdata2;
                
            case 5
                tempdata1                            = dataOut(:,:,3);
                tempdata2                            = dataOut(:,:,2);
                tempdata1(assignedNuclei==counterN)  = 255;
                tempdata2(assignedNuclei==counterN)  = 255;
                dataOut(:,:,3)                      = tempdata1;
                dataOut(:,:,2)                      = tempdata2;
                
            case 6
                tempdata1                            = dataOut(:,:,1);
                tempdata2                            = dataOut(:,:,2);
                tempdata3                            = dataOut(:,:,3);
                tempdata1(assignedNuclei==counterN)  = 255;
                tempdata2(assignedNuclei==counterN)  = 255;
                dataOut(:,:,1)                      = tempdata1;
                dataOut(:,:,2)                      = tempdata2;
                dataOut(:,:,3)                      = tempdata3;
                
        end
        %     % assigns tubules divided in 3
        %     currentLevel                        = 1+rem(counterN,3);
        %     tempdata                            = dataOut(:,:,currentLevel);
        %     tempdata(assignedNuclei==counterN)  =255;
        %     dataOut(:,:,currentLevel)           = tempdata;
    end
    
    %imagesc(dataOut)
    dataOut2=uint8(dataOut2);

