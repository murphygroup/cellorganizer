function status = tz_chkreg(img)
%TZ_CHKREG Check if regions in image have overlapped labels.
%   STATUS = TZ_CHKREG(IMG) returns 1 if there are pixels that have the
%   same gray level but are separated by background. Otherwise it returns 
%   0.
%   
%   See also

%   23-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

bwimg = img>0;
newLabeledImage = bwlabel(bwimg);

connectedObjectsNumber = max(newLabeledImage(:));
labeledNumber = max(img(:));

labelMatrix = zeros(connectedObjectsNumber,labeledNumber);

for i=1:connectedObjectsNumber
    labels = unique(img(find(newLabeledImage==i)));
    labelMatrix(i,labels) = 1;
end

status = all(sum(labelMatrix,1)<=1);