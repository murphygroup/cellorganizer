function cell_shape_model = train_cell_shape_model3( celltemppath,savepath )
% Train cell shape model using 3D Hela images

% Author: Devin Sullivan - adapted from Tao Peng's train_cell_shape_model2.m
%
% Copyright (C) 2011-2013 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
%Changes prior to 6/13/13 refactoring:
% July 26, 2012 I. Cao-Berg Added a statement where it returns an empty model
%               when the ratios of radii is empty
% May 15, 2013 I. Cao-Berg Updated method to support wildcards
%
%%
% June 13, 2013 D. Sullivan refactored to use the per-cell parameters
%                           precomputed using cellfit_percell.m
%%
%Changes post 6/13/13 refactoring:
%
%
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
% savepath = [pwd filesep 'temp' filesep 'cell_shape_eigen'];

warning('This function is depricated. Please see param2cell_ratio_model.m')

%if no temp path was specified, try default of pwd/temp/cell_shape_eigen
if isempty(celltemppath)
    celltemppath = [pwd filesep 'temp' filesep 'cell_shape_eigen'];
end
%next check if the temp dir exists
if ~exist(celltemppath,'dir')
    error('Unable to find per-cell models, please check celltempdir.');
end

cell_shape_model_file = [savepath filesep 'cell_shape_model.mat'];

if ~exist(cell_shape_model_file, 'file')

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
    [rad_ratio,cellnucratios,nucbottomslices] = extract_radius_ratio_old( celltemppath );

    %icaoberg 26/7/2012
    if isempty( rad_ratio )
       cell_shape_model = [];
       return
    end

    %rad_ratio(102,:) = [];
    %icaoberg 26/11/2013
    rad_ratio( (rad_ratio>1) ) = 1;

    keepinds = ~any(isnan(rad_ratio),2);
    
    if any(~keepinds)
        warning([num2str(sum(~keepinds)) ' cell shapes were removed due to improper segmentation.']);
    end
    % rad_ratio(any(isinf(rad_ratio),2),:) = [];

    %  rad_ratio(isinf(rad_ratio)) = 1;
    % tic
    if ~exist([savepath filesep 'pca_result.mat'],'file')
        [mu,coeff,score,latent] = eff_PCA(rad_ratio(keepinds,:));
        save([savepath filesep 'pca_result.mat'],'mu','coeff','score','latent')
    else
        load([savepath filesep 'pca_result.mat']);
    end
    % toc


    %% Cell shape model
    cell_shape_model = struct('name',[],'meanShape',[],'modeShape',[],'eigens',[]);
    cell_shape_model.name = 'RREigen';
    cell_shape_model.meanShape = struct('const',mu);

    try
        cell_shape_model.modeShape = struct('const',coeff(:,1:50));
        cell_shape_model.eigens = struct('stat',latent(1:50));
    catch err
        cell_shape_model.modeShape = struct('const',coeff(:,1:end));
        cell_shape_model.eigens = struct('stat',latent(1:end));
    end

    cell_shape_model.cellnucratio.stat = ml_estpdf(cellnucratios,struct('name','norm'));
    cell_shape_model.nucbottom.stat = ml_estpdf(nucbottomslices,struct('name','norm'));

    save(cell_shape_model_file,'cell_shape_model')

else
    disp(['Cell shape model file found. Loading ' cell_shape_model_file]);
    load(cell_shape_model_file)
end
