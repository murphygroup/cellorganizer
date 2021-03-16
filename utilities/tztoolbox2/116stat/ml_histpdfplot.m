function h = tz_histpdfplot(x,f,param)
%ML_HISTPDFPLOT Plot the histogram along with a pdf
%   H = ML_HISTPDFPLOT(X,F) shows the histogram of the vector X and the
%   [pdf] F. It also returns the handle of the plot.
%   
%   H = ML_HISTPDFPLOT(X,F,PARAM) specifies how to plot the figure by
%   PARAM: 
%       'nbin' - number of bins of the histogram (default 30)
%       'EdgeColor' - edge color of the bars (default [1 1 1])
%       'FaceColor' - face color of the bars(default [0.5 0.5 0.5])
%       'tz_plotfun' - parameters for TZ_PLOTFUN. The default plot parameters
%       are {'k-','lineWidth',2}.
%
%   Notice: if F is empty, only histogram will be plotted.
%
%   See also

%   20-Jul-2006 Initial write T. Zhao
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

if nargin < 1
    error('1 or 2 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

if isfield(f,'transform')
    x = ml_evalfun(x,f.transform);
end

param = ml_initparam(param,struct('nbin',30,'EdgeColor',[0 0 0], ...
    'FaceColor',[0.5,0.5,0.5],'tz_plotfun',struct( ...
    'range',[min(x),max(x)],'title','', ...
    'plot',{{'k-','lineWidth',2}})));

[n,y] = hist(x,param.nbin);
h = bar(y,n/length(x)/(y(2)-y(1)),'hist');
set(h,'EdgeColor',param.EdgeColor);
set(h,'FaceColor',param.FaceColor);

if ~isempty(f)
    hold on
    if isfield(f,'transform')
        f = rmfield(f,'transform');
    end
    tz_plotfun(ml_pdf2fun(f),param.tz_plotfun);
    hold off
end
