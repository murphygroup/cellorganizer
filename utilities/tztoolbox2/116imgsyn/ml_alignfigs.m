function h = ml_alignfigs(figs,param)
%ML_ALIGNFIGS
%   H = ML_ALIGNFIGS(FIGS)
%   
%   H = ML_ALIGNFIGS(FIGS,PARAM)
%   
%   See also

%   03-Apr-2007 Initial write T. Zhao
%   Copyright (c) 2007 Murphy Lab
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


if nargin < 1
    error('1 or 2 arguments are required');
end

param = ml_initparam(param,struct('size',[],'overwrite','no', ...
        'ncol',3,'nrow',1,'position',[.05 .05], ...   
        'marginal',0.1,'colspc',0.05,'rowspc',0.05, ...
        'text',{{'units','pixels','FontSize',18, ...
        'FontWeight','bold'}}, ...
        'savedir',[],'align',1, ...
        'axes',{{'FontSize',12,'LineWidth',2, ...
        'FontWeight','bold','units','pixels'}}, ...
        'xlabels',{{}},'ylabels',{{}}));

if isempty(param.size)
    param.size = [300 300];
end



figaxes = {};
for i=1:length(figs)
    open(figs{i});
    figaxes{i} = get(gca,'Children');
end

h = figure,
%axes('position',[0.1 0.1 1 1]);
set(h,'position',[80 50 param.size]);

width = (1-2*param.marginal-(param.ncol-1)*param.colspc)/param.ncol;
height = (1-2*param.marginal-(param.nrow-1)*param.rowspc)/param.nrow;
switch param.align
    case 2  
        offset = [param.marginal,1-param.marginal-0.22];
    case 1
        offset = [param.marginal,param.marginal];
    otherwise
        offset = [param.marginal,1-param.marginal-param.align];
end

axsize = [width height].*param.size;

for k=1:2
    if param.position(k)>0 & param.position(k)<1
        param.position(k) = axsize(k)*param.position(k);
    else
        param.position(k) = param.position(k);
    end
end

for i=1:length(figaxes)
    subplot('position',[offset width height]);
    copyobj(figaxes{i},gca);
    box
    if ~isempty(param.axes)
        set(gca,param.axes{:});
    end
    if ~isempty(param.xlabels)
        if ~isempty(param.xlabels{i})
            xlabel(param.xlabels{i});
        end
    end
    if ~isempty(param.ylabels)
        if ~isempty(param.ylabels{i})
            ylabel(param.ylabels{i});
        end
    end

    %set(gca,'FontSize',12);
    %set(gca,'LineWidth',2);
    %set(gca,'FontWeight','bold');
    %set(gca,'units','pixels');

    if ~isempty(param.label)
        text('String',['(' param.label ')'], ...
             'position',[param.position(1),axsize(2)-param.position(2)], ...
        param.text{:});
        param.label = char(param.label+1);
    end

    if mod(i,param.ncol)==0
        offset(1) = param.marginal;
        offset(2) = offset(2)-height-param.rowspc;
    else
        offset(1) = offset(1)+width+param.colspc;
    end
end


