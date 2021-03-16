function [img2,overlap] = tz_imaddobj(object,img,dist,center,option)
%TZ_IMADDOBJ Add an object to an image according to distance.
%   IMG2 = TZ_IMADDOBJ(OBJECT,IMG,DIST) adds the object OBJECT into the
%   image IMG according to the distance to the center of the image. OBJECT
%   has three column with the form [x,y,gray level].  It is the same as
%   IMG2 = TZ_IMADDOBJ(OBJECT,IMG,DIST,[])
%   
%   IMG2 = TZ_IMADDOBJ(OBJECT,IMG,DIST,CENTER) customizes the center to
%   the 2-element vector CENTER instead of the center of the image.
%
%   IMG2 = TZ_IMADDOBJ(OBJECT,IMG,DIST,CENTER,OPTION) specifies some 
%   special ways of patching objects by OPTION, which is a cell array of
%   a series of parameters:
%   ------------------------------------------------------------------
%   parameter      | data type |  description         | Default value
%   ------------------------------------------------------------------
%    overlap*1     | numeric   | overlapping or not   |      0
%   position*2     | string    | position pick        |   'random'   
%   object scale*3 | string    | adjust object values |   'original'
%   ------------------------------------------------------------------
%   *1 It has two values: 0 or 1. If it is 1, the position of object can
%   be automatically adjusted to avoid overlapping. If it is 2, the 
%   objects are not allowed to connect with each other.
%   
%   *2 It has two values:
%       'random' - randomly pick a position
%       'optimal' - pick a position farthest away from existed objects.
%   *3 It has two values:
%       'original' - orginal gray levels of objects after adding to the
%           image.
%       'normal' - normalized gray levels of objects
%
%   If a parameter in OPTION is emtpy or nonexistent, the default value
%   will be used.
%
%   IMG2 = TZ_IMADDOBJ(OBJECT,IMG,DIST,CENTER,{1}) can adjust the position
%   of the object automaticallly to avoid overlap.
%
%   [IMG2,OVERLAP] = TZ_IMADDOBJ(...) will also return the overlapped 
%   number of pixels between the new object and the old objects in the
%   image.

%   13-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('At least 3 arguments are required')
end

if nargin < 4
    center = [];
end

if nargin <5 
    nooverlap = [];
    positionPick = [];
    objectsc = [];
else
    nInput = length(option);
    for i=nInput+1:3
        option{i} = [];
    end

    nooverlap = option{1};
    positionPick = option{2};
    objectsc = option{3};
end

if isempty(nooverlap)
    nooverlap = 0;
end

if isempty(positionPick)
    positionPick = 'optimal';
end

if isempty(objectsc)
    objectsc = 'original';
end

if isempty(center)
    center = round(size(img)/2);
end

centralImg = zeros(size(img));
centralImg(center(1),center(2)) = 1;
centralMap = bwdist(centralImg);

validInd = find(abs(centralMap - dist)<=2);

switch positionPick
    case 'optimal'
        objectMap = bwdist(img>0);
        objectDists = objectMap(validInd);
        [sortedDist,distRank] = sort(objectDists,'descend');
    case 'random'
        distRank = randperm(length(validInd));
end

bestInd = distRank(1);
% [maxDist,bestInd] = max(objectDists);

[bestPosX,bestPosY] = ind2sub(size(img),validInd(bestInd));
object = ml_addrow(object,-[tz_calcobjcof(object),0]);
singleObjImg = tz_obj2img(ml_addrow(object,[bestPosX,bestPosY,0]), ...
    size(img),{'2d','og'});


if nargout==2 | nooverlap>0
    overlap = sum( sum( (singleObjImg>0) & (img>0) ) );
    maxdist = 100; 
    curdist = 2;
    olddist = curdist;
    while nooverlap == 1 & overlap>0 & curdist < maxdist
        distRank(1) = [];
%         objectDists(bestInd) = [];
        
        if isempty(distRank)
%             [maxDist,bestInd] = max(objectDists);
            
%         else
            curdist = curdist+2
            validInd = find(abs(centralMap - dist)<=curdist & ...
                abs(centralMap - dist)>olddist);
            
            switch positionPick
                case 'optimal'
                    objectDists = objectMap(validInd);
                    [sortedDist,distRank] = sort(objectDists,'descend');
                case 'random'
                    distRank = randperm(length(validInd));
            end
            olddist = curdist;
        end
        bestInd = distRank(1);
        [bestPosX,bestPosY] = ind2sub(size(img),validInd(bestInd));
        object = ml_addrow(object,-[tz_calcobjcof(object),0]);
        singleObjImg = tz_obj2img(ml_addrow(object,[bestPosX,bestPosY,0]), ...
            size(img),{'2d','og'}); 
        switch nooverlap
        case 1
            overlap = sum( sum( (singleObjImg>0) & (img>0) ) );
        case 2
            overlap = sum( sum( (filter2(ones(3,3), ...
                singleObjImg,'same')>0) & (img>0) ) );
        end
    end
end

if overlap>0 & nooverlap>0
    img2 = img;
else
    switch objectsc
    case 'original'
        img2 = singleObjImg+img;
    case 'normal'
        img2 = mat2gray(singleObjImg)+mat2gray(img);
    end
end
