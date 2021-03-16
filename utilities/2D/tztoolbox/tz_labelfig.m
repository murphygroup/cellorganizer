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

% Copyright (C) 2007  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

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