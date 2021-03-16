function model = train_cell_shape_model2( dnaImagesDirectoryPath, ...
    cellImagesDirectoryPath, ...
    proteinImagesDirectoryPath, ...
    param )
% Train cell shape model using 3D Hela images

% Author: Tao Peng
%
% Copyright (C) 2011-2013 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% July 26, 2012 I. Cao-Berg Added a statement where it returns an empty model
%               when the ratios of radii is empty
% May 15, 2013 I. Cao-Berg Updated method to support wildcards
%
%%
% June 10 2013 D. Sullivan Changed from extract_radius_ratio.m which now
%                          precomputes radius ratio statistics per-cell
%                          while extract_radius_ratiomodel builds the whole
%                          pattern model
%%
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

%%  Processing images and extract cell boundaries as cell codes
% Note this step can be parallelized
savepath = [pwd filesep 'temp' filesep 'cell_shape_eigen'];
if ~exist(savepath,'dir')
    mkdir(savepath)
end

% %icaoberg 5/15/2013
% try
%     maskImagesDirectoryPath = param.masks;
% catch
%     maskImagesDirectoryPath = '';
% end
% find_cell_codes( dnaImagesDirectoryPath, ...
%     cellImagesDirectoryPath, ...
%     proteinImagesDirectoryPath, ...
%     maskImagesDirectoryPath, savepath, param );

%% Pack up results and run PCA analysis


% June 10 2013 D. Sullivan Changed from extract_radius_ratio.m which now
%                          precomputes radius ratio statistics per-cell
%                          while extract_radius_ratiomodel builds the whole
%                          pattern model
[rad_ratio,cellnucratios,nucbottomslices] = extract_radius_ratiomodel( savepath );

%icaoberg 26/7/2012
if isempty( rad_ratio )
   model = [];
   return
end

%rad_ratio(102,:) = [];
rad_ratio(rad_ratio>1) = 1;

[mu,coeff,score,latent] = eff_PCA(rad_ratio);
save([savepath '/pca_result.mat'],'mu','coeff','score','latent')

%% Cell shape model
model = struct('name',[],'meanShape',[],'modeShape',[],'eigens',[]);
model.name = 'RREigen';
model.meanShape = struct('const',mu);

try
    model.modeShape = struct('const',coeff(:,1:50));
    model.eigens = struct('stat',latent(1:50));
catch err
    model.modeShape = struct('const',coeff(:,1:end));
    model.eigens = struct('stat',latent(1:end));
end

model.cellnucratio.stat = ml_estpdf(cellnucratios,struct('name','norm'));
model.nucbottom.stat = ml_estpdf(nucbottomslices,struct('name','norm'));

save([savepath '/cell_shape_model.mat'],'model')
