function [cellTubules_L,cellBody_L_Complete]=allocateTubules(cellBody_L,cellProtrusions,cellTubules,cellNoNuclei)

%% join cellBody with cellProtrusions
% Some protrusions will be close to just one cell, allocate those first
[cellProtrusions_L,numProt]             = bwlabel(cellProtrusions);
cellBody_LD                             = imdilate(cellBody_L,ones(5));

numCells                                = max(cellBody_L(:));
%%
cellProtrusions_L2 = zeros(size(cellProtrusions));
for counterProt = 1:numProt
    currentProt             = (cellProtrusions_L==counterProt);
    currentOverlap          = cellBody_LD.*currentProt;
    classProt               = unique(currentOverlap);
    numClassesProt          = numel(classProt);
    if numClassesProt==2
        % only two elements, [0 X] allocate to the second element
        cellProtrusions_L2  = cellProtrusions_L2 + classProt(2)*currentProt;
    else
        % more than one element, should allocate to whichever is closest
        % ... or has more of it 
        numElements = zeros(numClassesProt-1,1);
        for counterClasses = 1:numClassesProt-1
            numElements(counterClasses) = sum(currentOverlap(:)==classProt(counterClasses+1));
        end
        [~,indMax]= max(numElements);
        cellProtrusions_L2  = cellProtrusions_L2 + classProt(indMax+1)*currentProt;
        
    end
end
%%

%%
cellBody_L_Complete     = cellProtrusions_L2+cellBody_L;
cellBody_L_CompleteD    = imdilate(cellBody_L_Complete,ones(35));

%%
%imagesc(cellBody_L_CompleteD.*cellTubules)

%% tubule allocation
% There are several ways in which the tubules can be allocated:
% 0) allocate a label to each tubule separately
% 1) allocate to the one that belongs the most, i.e. closest to the cell
% 2) allocate a new class of the sum of the classes when overlap
% 3) allocate as cell/contact/no cell

[cellTubules_L0,numTubules]          = bwlabel(cellTubules.*(cellBody_L_Complete==0));
cellTubules_L1                      = zeros(size(cellTubules_L0));
cellTubules_L2                      = zeros(size(cellTubules_L0));
cellTubules_L3                      = zeros(size(cellTubules_L0));

for counterTub = 1:numTubules
    currentTub                     = (cellTubules_L0==counterTub);
    currentOverlap                  = cellBody_L_CompleteD.*currentTub;
    classTub                       = unique(currentOverlap);
    if numel(classTub)==2
        % only two elements, [0 X] allocate to the second element
        cellTubules_L1          = cellTubules_L1 + classTub(2)*currentTub;
        cellTubules_L2          = cellTubules_L2 + classTub(2)*currentTub;
        cellTubules_L3          = cellTubules_L3 + 1*currentTub;        
    elseif numel(classTub)>2
        % more than one element, should allocate to whichever is closest
        % ... or has more of it 
        for counterClasses = 1:numel(classTub)-1
            numElements(counterClasses) = sum(currentOverlap(:)==classTub(counterClasses+1));
        end
        [~,indMax]= max(numElements);
        % allocate to the one that belongs the most, i.e. closest to the
        % cell
        cellTubules_L1  = cellTubules_L1 + classTub(indMax+1)*currentTub;
        % or allocate a new class of the sum of the classes
        cellTubules_L2  = cellTubules_L2 + sum(classTub(2:end))*currentTub;
        cellTubules_L3  = cellTubules_L3 + 2*currentTub;
    else
        % tubule does not overlap with a cell, leave at 1 for the time being
         cellTubules_L1  = cellTubules_L1 +(numCells+2)*currentTub;
         cellTubules_L2  = cellTubules_L1 +(numCells+2)*currentTub;
         cellTubules_L3  = cellTubules_L3 + 3*currentTub;
    end
end

%% remove overlap with cells and save in a single variable
%imagesc(cellBody_L_Complete+cellTubules_L2)

cellTubules_L = {cellTubules_L0,cellTubules_L1,cellTubules_L2,cellTubules_L3};

%%
%qq=1;