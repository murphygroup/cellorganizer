function img2 = ml_imaddtext(img,txt,param)
%ML_IMADDTEXT Add text to an image.
%   IMG2 = ML_IMADDTEXT(IMG,TXT) returns an image that contains the text
%   TXT. 
%   
%   IMG2 = ML_IMADDTEXT(IMG,TXT,PARAM) customizes how to add the text by
%   specifying the structure PARAM with the following fields:
%       'font' - font of the text. Default 'courier'.
%       'fontsize' - size of the text. Default 18.
%       'color' - color of the text. Currently only two colors are
%           supported:
%               'w' : white (default)
%               'k' : black
%           It does not work for a constant image.
%       'position' - the position of the text. It is relative to the
%           left-top corner of IMG. first element corresponds to row and
%           the second element corresponds to column. If the value of an
%           element is greater than 1, then it is the number of pixels of
%           the offset of the text from the image corner. If the value of
%           an element is between 0 and 1, it is the ratio of the image
%           dimension to move the text. Default [.05 .05].
%
%   Notice: This function requires the graphical interface to work.
%
%   See also

%   26-Oct-2006 Initial write T. Zhao
%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 2
    error('2 or 3 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('fontsize',18,'color','w', ...
    'position',[.05 .05],'font','courier'));

h = figure('color','white','units','normalized','position',[.3 .3 .3 .3]);
set(gca,'units','pixels','position',[5 5 200 200], ...
    'visible','off')
text('units','pixels','position',[30 30], ...
    'fontsize',param.fontsize,'fontweight','bold','string',txt,'Color','k');
axis off

textimg = getframe(gca);
close(h)

% Extract the cdata
textimg = textimg.cdata(:,:,1)==0;
textimg = ml_imcropbg(textimg);
textobj = tz_img2obj(textimg);

if strcmp(param.color,'k')
    textobj(:,3) = min(img(:));
end
if strcmp(param.color,'w')
    textobj(:,3) = max(img(:));
end

for k=1:2
    if param.position(k)>0 & param.position(k)<1
        textpos(k) = size(img,k)*param.position(k);
    else
        textpos(k) = param.position(k);
    end
end

img2 = ml_imaddobj2(img,textobj, ...
    struct('pos',textpos,'posref','corner'));
