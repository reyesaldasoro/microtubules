function  [assignedTubules,finalCells,distFromCell]         = assignTubulesToCells(microTubules_L,nuclei_L,cellBody_L)

% find dimensions of the input data
[rows,cols]                                 = size(microTubules_L);
numMicroTubules                             = max(microTubules_L(:));
numNuclei                                   = max(nuclei_L(:));

    %% Assign tubules to a single cell based on the nuclei but distance TO CELL BODY
    
    %regionsCells                           = watershed(-nuclei);
    %[microTubules_L,numMicroTubules]        = bwlabel(microTubules);
    microTubules                        = (microTubules_L>0);
        
    
    nucleiDist(rows,cols,numNuclei+1)       = 0;
    distFromCell(rows,cols,numNuclei)       = 0;
    finalCells(rows,cols,numNuclei)         = 0;
    %%
    for counterN=1:numNuclei
        % The distance should be calculated from the cell Body, but cell bodies can
        % merge so split with nuclei if necessary
        % counterN=1;
        currentNuclei                       = (nuclei_L==counterN);
        currentCellIn                       = unique(cellBody_L(currentNuclei));
        currentCellIn(currentCellIn==0)     = [];
        %try
        currentCell                         = (cellBody_L==currentCellIn);
        %catch
        %    qqq=1;
        %end
        testNumbNuclei                      = unique(nuclei_L(currentCell));
        
        if numel(testNumbNuclei)==2
            finalCells(:,:,counterN)         = currentCell;
            distFromCell(:,:,counterN)      = bwdist(currentCell);
            nucleiDist(:,:,counterN)        = distFromCell(:,:,counterN).*(microTubules);
        else
            %split the cell with the nuclei
            cellSplit                       = watershed(-ismember(nuclei_L,testNumbNuclei(2:end)));
            cellSplitToKeep                 = unique(cellSplit(currentNuclei));
            finalCells(:,:,counterN)        = currentCell.*(cellSplit==cellSplitToKeep);
            distFromCell(:,:,counterN)      = bwdist(finalCells(:,:,counterN));
            nucleiDist(:,:,counterN)        = distFromCell(:,:,counterN).*microTubules;
        end
        
    end
    
    nucleiDist(:,:,numNuclei+1)             = -(1-microTubules)+1000*microTubules;
    [maxDistT,correspNuclei]                = min(nucleiDist,[],3);
    %%
    % kkk=4;
    % imagesc(nucleiDist(:,:,kkk)+100*(cellBody)+100*(nuclei_L==kkk))
    
    %%
    assignedTubules(rows,cols)               = 0;
    for counterMT =1:numMicroTubules
        %assign every microtubule to a cell body, do this by checking the median distance
        %of the microtubule and then assigning to a cell, BUT discard those that are far
        %from cells as Cells without nucleus have been discarded
        if mean(maxDistT(microTubules_L==counterMT))<50
        assignedTubules(microTubules_L==counterMT) =median(correspNuclei(microTubules_L==counterMT));
        end
    end

       