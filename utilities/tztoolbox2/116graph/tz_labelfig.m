function tz_labelfig(h,xl,yl,tt)
%TZ_LABELFIG Label a figure and make it suitable for publication.
%   TZ_LABELFIG(H) make the figure with handle H suitable for publication.
%   
%   TZ_LABELFIG(H,XL) labels the x axis by string XL.
%   
%   TZ_LABELFIG(H,XL,YL) labels the y axis by string YL.
%   
%   TZ_LABELFIG(H,XL,YL,TT) adds a title TT into the figure.
%
%   See also

%   24-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('At least 1 arguments are required')
end

p=get(h);
set(p.Children,'FontSize',12);
set(p.Children,'LineWidth',2);
set(p.Children,'FontWeight','bold');

if exist('xl','var')
    xlabel(['\bf\fontsize{14} ' xl]);
end

if exist('yl','var')
    ylabel(['\bf\fontsize{14} ' yl]);
end

if exist('tt','var')
    title(['\bf\fontsize{14} ' tt]);
end