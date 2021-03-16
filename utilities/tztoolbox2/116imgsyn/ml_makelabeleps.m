function savefile = ml_makelabeleps(fig,param)
%ML_MAKELABELEPS Make a labled eps file from a figure.
%   SAVEFILE = ML_MAKELABELEPS(FIG,PARAM) put a label in the figure FIG, which
%   is a handle of a figure or the file path of figure. The eps file will be
%   save in SAVEFILE, which is determeined by the fields of the stucture PARAM.
%   See ML_MAKELABELFIG for more details. Other fields of PARAM are:
%       'overwrite' - see ML_MAKELABELFIG.
%       'imgsize' - size of the eps image. Default [200 200].
%       'text' - a cell array of parameters for the TEXT function.
%   
%   See also

%   30-Oct-2006 Initial write T. Zhao
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
    error('Exactly 2 arguments are required');
end

param = ml_initparam(param,struct('imgsize',[],'overwrite','no', ...
        'position',[.05 .05], ...                          
        'text',{{'units','pixels','FontSize',18, ...
        'FontWeight','bold'}}, ...
        'savedir',[]));

separator = '.';

if ~isempty(param.savedir)
    savefile = [param.savedir filesep param.paperid separator ...
            num2str(param.figidx) separator param.label '.eps'];
else
    savefile = [];
end

if exist(savefile,'file')
    if strcmp(param.overwrite,'yes')
        warning(['The file ' savefile ' is overwritten']);
    else
        warning(['The file ' savefile ' already exists. Skip this ' ...
                            'file']);
        return;
    end
end

if isempty(param.imgsize)
    param.imgsize = [300 300];
end

for k=1:2
    if param.position(k)>0 & param.position(k)<1
        param.position(k) = param.imgsize(k)*param.position(k);
    else
        param.position(k) = param.position(k);
    end
end

if ischar(fig)
    fig = open(fig);
end
set(gca,'units','pixels','position',[80 50 param.imgsize]);

if ~isempty(param.label)
    text('String',['(' param.label ')'], ...
         'position',[param.position(1),param.imgsize(2)-param.position(2)], ...
         param.text{:});
end

if ~isempty(savefile)
    print(gcf,'-deps',savefile);
end

