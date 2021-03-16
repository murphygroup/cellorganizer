function h = ml_normplot(y,param)
%ML_NORMPLOT Plot the histogram along with the normal distribution.
%   H = ML_NORMPLOT(Y) shows the histogram of the vector Y and the
%   esitamited normal distribution of Y. It also returns the handle of the
%   bar plot.
%   
%   H = ML_NORMPLOT(Y,PARAM) specifies how to plot the figure by PARAM:
%       'nbin' - number of bins of the histogram (default 30)
%       'EdgeColor' - edge color of the bars (default [1 1 1])
%       'FaceColor' - face color of the bars(default [0.5 0.5 0.5])
%       'plot' - a cell array secifying plot options (default
%           {'k-','lineWidth',2})
%       'showpdf' - 'yes' to show the estimated pdf (default). 'no' for not 
%           showing. 
%   
%   See also

%   28-Jun-2006 Initial write T. Zhao
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


if nargin < 1
    error('1 or 2 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('nbin',30,'EdgeColor',[1 1 1]-1, ...
    'FaceColor',[0.5,0.5,0.5],'plot',{{'k-','lineWidth',2}},'showpdf','yes'));

[n,x] = hist(y,param.nbin);
h = bar(x,n/length(y)/(x(2)-x(1)),'hist');
set(h,'EdgeColor',param.EdgeColor);
set(h,'FaceColor',param.FaceColor);

if strcmp(param.showpdf,'yes')
    mu = mean(y);
    sigma = std(y);
    p = normpdf(x,mu,sigma);
    hold on
    plot(x,p,param.plot{:});
    hold off
end
