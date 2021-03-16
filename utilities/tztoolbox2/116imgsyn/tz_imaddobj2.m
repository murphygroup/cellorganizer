function img2 = tz_imaddobj2(img,obj,param)
%TZ_IMADDOBJ2 Add an object into an image. Obsolete. See ML_IMADDOBJ2.
%   IMG2 = TZ_IMADDOBJ2(IMG,OBJ) adds an [object] into the [image] IMG and
%   returns the new image. The position of OBJ is determined by its 
%   coordinates.
%   
%   IMG2 = TZ_IMADDOBJ2(IMG,OBJ,PARAM) customizes different ways of putting
%   OBJ. PARAM can have the following fields:
%   ------------------------------------------------------------------
%        name      | data type |  description           | Default value
%   ------------------------------------------------------------------
%      pos         | [point]   | position of OBJ        |      []
%     method       | string    | how to add OBJ         |   'replace' 
%      posref      | string    | how to specify object  |   'cof'
%                  |           | position               |
%   ------------------------------------------------------------------
%
%   More details:
%       'pos' - if it is empty, the position of the oject will be
%           determined by its coordinates.
%       'method' - it has two values. 'replace' means that the the values
%           in IMG will be replaced by OBJ intensities. 'add' means that
%           the OBJ intensities will be add into IMG.
%       'posref' - it is only applicable when 'pos' is not empty. If it is
%           'cof', then 'pos' will be where the COF of OBJ is located in
%           the image IMG. If it is 'corner', then 'pos' will be where the
%           left top coner of the bounding box of OBJ is located.
%   
%   See also

%   26-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_imaddobj2','ml_imaddobj2'));

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param, ...
    struct('method','replace','pos',[],'posref','cof'));

if ~isempty(param.pos)
    switch param.posref
        case 'cof'
            offset = tz_calcobjcof(obj);
        case 'corner'
            box = ml_boundbox(obj(:,1:2));
            offset = box(1,:);
    end

    obj(:,1:2) = round(ml_addrow(obj(:,1:2),-offset+param.pos));
end

img2 = img;
imageSize = size(img2);
if length(imageSize)>2
    imageSize = imageSize(1:2);
end

idx = ml_objinimg(obj,imageSize);
if ~isempty(idx)
    warning('The object is out of the image range!')
    obj(idx,:) = [];
end

objidx = sub2ind(imageSize,obj(:,1),obj(:,2));

if ~isempty(obj)
    switch param.method
        case 'replace'
            for k=1:size(img2,3)
                img3 = img2(:,:,k);
                img3(objidx) = obj(:,3);
                img2(:,:,k) = img3;
            end
            %img2(objidx) = obj(:,3);
        case 'add'
            for k=1:size(img2,3)
                img3 = img2(:,:,k);
                img3(objidx) = img3(objidx)+obj(:,3);
                img2(:,:,k) = img3;
            end
            %img2(objidx) = img2(objidx)+obj(:,3);
        otherwise
            error(['Unrecognized object adding method: ' param.method]);
    end
end
