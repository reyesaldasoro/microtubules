function []=analyseCellConditions(cellBody,cellNuclei,cellProtrusions)
%%

% Determine cells with 2 nuclei (clumps)
nucleiPresent_centroids                 =(bwmorph(cellNuclei,'shrink','inf'));
nucleiPresent_location                  = find(nucleiPresent_centroids);
nucleiPresent                           = cellBody(nucleiPresent_location);



