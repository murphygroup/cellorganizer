function model = train_cell_shape_model( imgpath, param )
% Train cell shape model using 3D Hela images

% Author: Tao Peng
%
% Copyright (C) 2011-2012 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% March 23, 2012 I. Cao-Berg Modified save path directory
% March 28, 2012 R.F. Murphy Change to 3D cell/nuclear ratios
% Jan 1, 2013 I. Cao-Berg Updated method to use display according to
% verbose flag
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

%icaoberg 1/30/2013
if nargin == 1
   param = [];
end

try
  verbose = param.verbose;
catch
  param.verbose = true;
end

try
  debug = param.debug;
catch
  param.debug = false;
end

savepath = [ pwd filesep 'temp' filesep 'cell_shape_eigen'];
if ~exist(savepath,'dir')
    mkdir(savepath)
end

if verbose; disp( 'Finding cell ratios' ); end
rad_ratio = find_cell_ratios( imgpath, savepath );
rad_ratio(rad_ratio>1) = 1;

[x,y]=find(isinf(rad_ratio));
x=unique(x);
for i=1:1:length(x)
   if verbose; disp(['Ignoring feature vector ' num2str(i)]); end
   rad_ratio = [rad_ratio(1:x-1,:,:); rad_ratio(x+1:end,:,:)];  
end

[x,y]=find(isnan(rad_ratio));
x=unique(x);

for i=1:1:length(x)
   if verbose; disp(['Ignoring feature vector ' num2str(i)]); end
   rad_ratio = [rad_ratio(1:x-1,:,:); rad_ratio(x+1:end,:,:)];  
end

%icaoberg april 2, 2012
if verbose; disp( 'Principal component analysis' ); end
[mu,coeff,score,latent] = eff_PCA(rad_ratio);

save([savepath '/pca_result.mat'],'mu','coeff','score','latent')

%cell shape model
if verbose; disp( 'Creating structure...' ); end
model = struct('name',[],'meanShape',[],'modeShape',[],'eigens',[]);
model.name = 'RREigen';
model.meanShape = struct('const',mu);

try
    model.modeShape = struct('const',coeff(:,1:20));
    model.eigens = struct('stat',latent(1:20));
catch err
    model.modeShape = struct('const',coeff(:,1:end));
    model.eigens = struct('stat',latent(1:end));
end

save([savepath '/cell_shape_model.mat'],'model')
