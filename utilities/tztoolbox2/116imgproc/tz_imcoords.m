function pts = tz_imcoords(imgsize,scale,offset)
%TZ_IMCOORDS Obsolete. See ML_IMCOORDS.
%   PTS = TZ_IMCOORDS(IMGSIZE) returns a 2x(MxN) matrix if IMGSIZE is
%   [M,N]. The coordinates are obtained column by column.
%
%   PTS = TZ_IMCOORDS(IMGSIZE,SCALE) will rescale the coordinates by
%   1/SCALE.
%
%   PTS = TZ_IMCOORDS(IMGSIZE,SCALE,OFFSET) will move the coordinates by 
%   OFFSET after scaling.

%   11-Sep-2005 Initial write T. Zhao

error(tz_genmsg('of','tz_imcoords','ml_imcoords'));

if nargin < 1
    error('At least 1 argument is required')
end

if nargin < 2
    scale = 1;
end

if nargin < 3
    offset = [0 0];
end

if length(scale)==1
    scale = zeros(1,2)+scale;
end

if(length(imgsize) ~= 2)
    error('The first argument must have exactly 2 elements');
end

xCoords = (1:imgsize(1))';
yCoords = 1:imgsize(2);

xPanel = xCoords( :,ones(imgsize(2),1) );
yPanel = yCoords(ones(imgsize(1),1),:);

pts = [xPanel(:)'/scale(1)+offset(1);yPanel(:)'/scale(2)+offset(2)];