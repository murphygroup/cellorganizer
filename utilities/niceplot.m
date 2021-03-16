function H = niceplot(xval,yval,xlabelstr,ylabelstr,titlestr,plottype,plotstyle)
%NICEPLOT Helper function to plot something that is readable and savable in a new figure window.
%
%Inputs:
%
% xval = vector of x values 
% yval = vector of y values
% xlabelstr = string for name of the x variable
% ylabelstr = string for name of the y variable
% titlestr = optional string to label the figure with a title
% plottype = string of plot type (scatter,plot,semilogx,semilogy,loglog)
% plotstyle = string e.g. '.' or '-B'
%
%Outputs:
% h = figure handle
%
%Author: Devin Sullivan 6/17/13

% Copyright (C) 2007-2013  Murphy Lab
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


if nargin<7
    plotstyle = '.';
elseif nargin==5
    plotstyle = '.';
    plottype = 'scatter';
elseif nargin==4
    plotstyle = '.';
    plottype = 'scatter';
    titlestr = '';
end


switch plottype

    case 'scatter'
        figure,scatter(xval,yval,100,plotstyle);
    case 'plot'
        figure,plot(xval,yval,plotstyle,'Linewidth',2,'MarkerSize',20);
    case 'semilogx'
        figure,semilogx(xval,yval,plotstyle,'Linewidth',2,'MarkerSize',20);
    case 'semilogy'
        figure,semilogy(xval,yval,plotstyle,'Linewidth',2,'MarkerSize',20);
    case 'loglog'
        figure,loglog(xval,yval,plotstyle,'Linewidth',2,'MarkerSize',20);
    otherwise
        disp('Sorry, unrecognized plot type');
        return
end
title(titlestr,'FontSize',20,'Fontname','Ariel')
xlabel(xlabelstr,'FontSize',20,'Fontname','Ariel');
ylabel(ylabelstr,'FontSize',20,'Fontname','Ariel');