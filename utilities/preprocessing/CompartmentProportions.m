function [totfluor,compfluor,propfluor] = ...
    CompartmentProportions(rawimg,compartimg,compartlist,param)
%This function takes a raw image an precomputed masks and returns the
%proportion of fluorescence in each of the compartments defined by
%DefineCompartments.m(cell,nuc,membrane)
%
%Inputs: 
%rawimg = a raw single cell image, this could be preprocessed(thresholded
%etc) but should have real values (e.g. not boolean unless you just want
%area/compartment)
%compartimg = an indexed image where each index is a compartment space, 0
%             is assumed to be the extracellular space which is discarded
%compartlist = a cell array naming the various compartments
%              both compartimg and compartlist can be obtained using 
%              DefineCompartments.m
%
%Outputs: 
%propfluor = the proportion of fluorescence in each compartment

indexes = unique(compartimg);
if (length(indexes)-1)~=length(compartlist)
    error('Number of compartments detected does not match number of compartments listed');
end

%check if passed a image path or actual raw image 
if isnumeric(rawimg)
    %already have an image, do nothing
else
    %probably have image path, try to read image in
    rawimg = ml_readimage(rawimg);
end

% %Only take the segmented part of the image 
% rawimg = rawimg(:,:,param.botslice:param.topslice);

%check the size of rawimg and masks are the same. if not, resize masks to
%fit the size of rawimg. 
%first check z for efficiency
if size(rawimg,3)~=size(compartimg,3)
    compartimg = tp_stretch3d(compartimg,size(rawimg,3),'nearest');
end
if any(size(rawimg)~=size(compartimg))
    %again, assume equal resolution in x and y
    compartimg = imresize(compartimg,size(rawimg,1)/size(compartimg,1),'nearest');
end


%for each compartment in the list, find the amount of fluorescence
%first mask with the compartments
allcomp = rawimg.*(compartimg>0);
%next get total fluorescence in the compartments
totfluor = sum(allcomp(:));
for i = 1:length(compartlist)
    currcompraw = (compartimg==indexes(i+1)).*rawimg;
    compfluor(i) = sum(currcompraw(:));
    propfluor(i) = compfluor(i)/totfluor;
end

%Lastly, if the user has specified a pattern mask they would like to find
%the fluorescence of, find it
%(this will no longer sum to 1)
%additional masks should be a cell array of masks the user wants info about
if isfield(param,'additionalmasks') && iscell(param.additionalmasks)
    %add counter of actual mask
    j = 1;
    for i = length(compartlist)+1:(length(param.additionalmasks)+length(compartlist))
        compfluor(i) = sum(param.additionalmasks{j}(:));
        propfluor(i) = compfluor(i)/totfluor;
        j = j+1;
    end
elseif isfield(param,'additionalmasks') && isnumeric(param.additionalmasks)
    i = i+1;
    compfluor(i) = sum(param.additionalmasks(:));
    propfluor(i) = compfluor(i)/totfluor;
end