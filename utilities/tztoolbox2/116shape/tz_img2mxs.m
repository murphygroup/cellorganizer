function [medaxis,width] = tz_img2mxs(img)
%TZ_IMG2MXS Convert an image into medial axis representation.
%   MEDAXIS = TZ_IMG2MXS(IMG) returns the medial axis of the image. This
%   function will take all pixels in IMG above intensity 0 as objects.
%   
%   [MEDAXIS,WIDTH] = TZ_IMG2MXS(...) also returns width.
%   
%   See also

%   01-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

[imgaxis,medaxis,width,borders]=tz_imaxis(img);