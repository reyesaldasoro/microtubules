function [cellTubules_L]=allocateTubules(cellBody_L,cellProtrusions,cellTubules)

%% join cellBody with cellProtrusions
% Some protrusions will be close to just one cell, allocate those first
[cellProtrusions_L,numProt]             = bwlabel(cellProtrusions);
cellBody_LD                             = imdilate(cellBody_L,ones(5));
%%
cellProtrusions_L2 = zeros(size(cellProtrusions));
for counterProt = 1:numProt
    currentProt             = (cellProtrusions_L==counterProt);
    currentOverlap          = cellBody_LD.*currentProt;
    classProt               = unique(currentOverlap);
    if numel(classProt)==2
        % only two elements, [0 X] allocate to the second element
        cellProtrusions_L2  = cellProtrusions_L2 + classProt(2)*currentProt;
    else
        % more than one element, should allocate to whichever is closest
        % ... or has more of it 
        for counterClasses = 1:numel(classProt)-1
            numElements(counterClasses) = sum(currentOverlap(:)==classProt(counterClasses+1));
        end
        [maxNum,indMax]= max(numElements);
        cellProtrusions_L2  = cellProtrusions_L2 + classProt(indMax+1)*currentProt;
        
    end
end
%%

imagesc(cellProtrusions_L2)
