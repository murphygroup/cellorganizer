function model = gmm_objsizefit_percell( objsize)
% Learn object size model for each location pattern

% Author: Tao Peng
% Edited: Ivan E. Cao-Berg
%
% Copyright (C) 2011-2016 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
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

% April 16, 2012 I. Cao-Berg Removed protein type from function call
% July 30, 2012 I. Cao-Berg Updated object size model to be a fit 
%                           from an exponential to a lognormal
% August 6, 2012 D. Sullivan major bug fix: added separate array for xbin 
%           values to prevent the x values from being overwritten before 
%           the z_x distribution is learned
% August 7, 2012 R.F. Murphy fix loop limit
%  
%
% Edited: D. Sullivan 8/6/2012 undid some of the previous changes such as
% binsizes. Go from 0-xobjsize==30 (arbitrarily chosen) and use a max_xbin
% = 15 (arbitrarily chosen).
%
% Edited: D. Sullivan 8/7/2012 reverted entirely to original version to fix
% issue with objects no longer appearing gaussian. Then increased range of
% x from 0-xobjsize==30 (arbitrarily chosen) and use a max_xbin
% = 15 (arbitrarily chosen).
%
% Edited: D. Sullivan 6/11/13 - refactored for per-cell computations. Also
%                             removed arbitrary max_xbin and max sizes for
%                             y and z in favor of 95% confidence. This is
%                             resolution invarient. 

 
% if ~isdir(savedir)
%     mkdir(savedir)
% end
% if ~exist([savedir filesep 'sizefits_' num2str(currfile) '.mat'],'file')
% 
% %D. Sullivan 6/11/13 removing resolution dependence, basing on 95%
% %   confidence interval
% %arbitrarily chosen to be the maximum size of xsigma to consider
% % max_xbin = 15;
% 
% file = ml_dir([datadir filesep '*' num2str(currfile) '.mat']);
% load([datadir filesep file{1}]);
objsize = real(objsize);

x = objsize(:,1);
y = objsize(:,2);
z = objsize(:,3);

%D. Sullivan 6/11/13 removing resolution dependence, basing on 95%
%   confidence interval X DIMENSION
tmpx = sort(x);
x95 = tmpx(1:floor(0.95*length(x)));
max_xbin = floor(x95(end));

if max_xbin == 0
    model = [];
    return
end

%D. Sullivan 6/11/13 removing resolution dependence, basing on 95%
%   confidence interval Y DIMENSION
ty = sort(y);
y95 = ty(1:floor(0.95*length(ty)));
optvars.alpha1y = y95(end);%(LAMP this = ~5.9)vs6.5
optvars.alpha2y = std(y95);%(LAMP this = ~1.3)vs1.4

%D. Sullivan 6/11/13 removing resolution dependence, basing on 95%
%   confidence interval Z DIMENSION
tz = sort(z);
z95 = tz(1:floor(0.95*length(tz)));
optvars.alpha1z = z95(end);
optvars.alpha2z = std(z95);

%icaoberg 7/30/2012
%lambda = expfit(x);
%model.x = struct('name','exp','beta',lambda);
%D. Sullivan, 6/11/13 now only fitting 95% confidence sizes
[parmhat,parmci] = lognfit( x95 );
model.x = struct('name', 'lognorm', 'mu', parmhat(1), 'sigma', parmhat(2) );

% Fit the y-x conditional normal distribution
Y_X = cell(max_xbin,1);
% actual_max_xbin = 1;
for i=1:max_xbin
    y_given_x = y(x>i-1&x<=i);
    Y_X{i} = y_given_x;
    %D. Sullivan 6/11/13 - added else in case there are no observed objects
    %in that bin
    if ~isempty(y_given_x)
        [mu(i),sigma(i)] = normfit(y_given_x);
%         actual_max_xbin = i;
    else
        mu(i) = 0;
        sigma(i) = 0;
    end
end

% dpsulliv 8/7/12 added max_xbin as a variable
% x = .5:14.5;
xbin = .5:max_xbin-0.5;

%D. Sullivan 6/11/13 removing resolution dependence, basing on 95%
%   confidence interval
% optvars.alpha1 = 6.5;
% optvars.alpha2 = 1.4;
% dpsulliv 8/7/12 added max_xbin as a variable
% u = mu(1:15);
u = mu(1:max_xbin);
xbin(u>=optvars.alpha1y) = [];
u(u>=optvars.alpha1y) = [];
optvars.beta1y = sum(-xbin.*log(1-u/optvars.alpha1y))/sum(xbin.^2);
% dpsulliv 8/7/12 added max_xbin as a variable
% u = sigma(1:15);
u = sigma(1:max_xbin);
% x = .5:14.5;
xbin = .5:max_xbin-0.5;
xbin(u>=optvars.alpha2y) = [];
u(u>=optvars.alpha2y) = [];
optvars.beta2y = sum(-xbin.*log(1-u/optvars.alpha2y))/sum(xbin.^2);

model.y_x = struct('name','norm',...
                    'a1',optvars.alpha1y,...
                    'b1',optvars.beta1y,...
                    'a2',optvars.alpha2y,...
                    'b2',optvars.beta2y);

% Fit the z-x conditional normal distribution
%clear mu sigma
% dpsulliv 8/7/12 expanded the range to 30 
% Z_X = cell(20,1);
% for i=1:20
Z_X = cell(max_xbin,1);
for i=1:max_xbin
    z_given_x = z(x>i-1&x<=i);
    Z_X{i} = z_given_x;
    if ~isempty(z_given_x)
        [mu(i),sigma(i)] = normfit(z_given_x);
    end
end

% dpsulliv 8/7/12 added max_xbin as a variable
% x = .5:14.5;
xbin = .5:max_xbin-0.5;
%D. Sullivan 6/11/13 removing resolution dependence, basing on 95%
%   confidence interval
% optvars.alpha1 = 3.5;
% optvars.alpha2 = 1.4;
% dpsulliv 8/7/12 added max_xbin as a variable
% u = mu(1:15);
u = mu(1:max_xbin);
xbin(u>=optvars.alpha1z) = [];
u(u>=optvars.alpha1z) = [];
optvars.beta1z = sum(-xbin.*log(1-u/optvars.alpha1z))/sum(xbin.^2);
% dpsulliv 8/7/12 added max_xbin as a variable
% u = sigma(1:15);
u = sigma(1:max_xbin);
% x = .5:14.5;
xbin = .5:max_xbin-0.5;
xbin(u>=optvars.alpha2z) = [];
u(u>=optvars.alpha2z) = [];
optvars.beta2z = sum(-xbin.*log(1-u/optvars.alpha2z))/sum(xbin.^2);

model.z_x = struct('name','norm',...
                    'a1',optvars.alpha1z,...
                    'b1',optvars.beta1z,...
                    'a2',optvars.alpha2z,...
                    'b2',optvars.beta2z);
