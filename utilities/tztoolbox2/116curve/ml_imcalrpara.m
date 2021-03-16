function [calpara,xdata,ydata,maxdist] = ml_imcalrpara(img,cut)
%ML_IMCALRPARA Calibrate the gray level of an object.
%   CALPARA = ML_IMCALRPARA(IMG) returns the calibration parameters of
%   the image [IMG] to make the object more uniform upon gray levels.
%   
%   CALPARA = ML_IMCALRPARA(IMG,CUT) lets users specify the cut level
%   for removing data, which means pixels with intensities lower than CUT
%   will be considered as background.
%   
%   [CALPARA,XDATA,YDATA,MAXDIST] = ML_IMCALRPARA(...) also returns data for
%   calibration. MAXDIST is the max distance of pixels inside the object to
%   object boundry.
%   
%   See also TZ_TEXCALRPARA

%   21-Aug-2006 Initial write T. Zhao
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

if ~exist('cut','var')
    cut=0.5;
end

distr=tz_getdistgray(img);
maxdist = max(distr(:,1));
xdata=(max(distr(:,1))-distr(:,1))/maxdist;
ydata=distr(:,2);

ydata(xdata<cut)=[];
xdata(xdata<cut)=[];

% [calpara,resid]=lsqcurvefit(@tz_projball,[10,1,max(xdata)], ...
%     xdata,ydata,[0,0,max(xdata)],[Inf,Inf,Inf]);
% [calpara,resid]=lsqcurvefit(@tz_projball,calpara, ...
%     xdata,ydata,[0,0,max(xdata)],[Inf,Inf,Inf]);

% [calpara,resid]=lsqcurvefit(@tz_projball,[10,1], ...
%     xdata,ydata,[0,0],[Inf,Inf]);
% [calpara,resid]=lsqcurvefit(@tz_projball,calpara, ...
%     xdata,ydata,[0,0],[Inf,Inf]);

calpara = ml_fitprojball(xdata,ydata);
%[calpara,resid]=lsqcurvefit(@tz_projball,calpara,xdata,ydata,[0,0,max(xdata)],[Inf,Inf,Inf]);