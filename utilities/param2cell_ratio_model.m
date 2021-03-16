function cell_shape_model = param2cell_ratio_model( cellparam , options)
% Train cell shape model using 3D Hela images

% Author: Devin Sullivan - adapted from Tao Peng's train_cell_shape_model2.m
%
% Copyright (C) 2011-2013 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% Adapted from train_cell_shape_model3.m - grj 10/13/15
%
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

nparam = length(cellparam);

for i = 1:nparam
    if ischar(cellparam{i})
        cellparam{i} = load(cellparam{i}, 'cell');
    end
end

if iscell(cellparam(1))
    cellparam = [cellparam{:}];
end

rad_ratio = cell(1, nparam);
cellnucratios = zeros(1, nparam);
nucbottomslices = zeros(1, nparam);

keepind = true(1, nparam);


for i = 1:nparam
    if isfield(cellparam(i).cell, 'rad_ratio')
        rad_ratio{i} = cellparam(i).cell.rad_ratio;
        cellnucratios(i) = cellparam(i).cell.cellnucheightratio;
        nucbottomslices(i) = cellparam(i).cell.nucbottomslice;
    else
        keepind(i) = false;
    end
end

rad_ratio = vertcat(rad_ratio{keepind});
cellnucratios = cellnucratios(keepind)';
nucbottomslices = nucbottomslices(keepind)';


%     [rad_ratio,cellnucratios,nucbottomslices] = extract_radius_ratio( celltemppath );

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
    
    [mu,coeff,score,latent] = eff_PCA(rad_ratio(keepinds,:));
%     if ~exist([savepath filesep 'pca_result.mat'],'file')
%         
%         save([savepath filesep 'pca_result.mat'],'mu','coeff','score','latent')
%     else
%         load([savepath filesep 'pca_result.mat']);
%     end
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
    
    cell_shape_model.type = 'ratio';
end

function [rad_ratio,cellnucratios,nucbottomslices] = extract_radius_ratio(cellcodepath)
% From the cellcodes, extract the radius ratio of nuclear radius over cell
% radius under cylindrical systerm
%
% Adapted from extract_radius_ratio_old

% Constants
n = 360;
H = -pi:(2*pi/360):pi;
k = 0;

% Calculate radius ratios and combining them
distratio = []; cellnucratios = []; nucbottomslices = [];
celllist = ml_dir([cellcodepath filesep 'cellcodes_*']);

for i = 1:length(celllist)
    fname = [cellcodepath filesep celllist{i}];
    if exist(fname,'file')
        k = k + 1;
        load(fname)
        baseplane(k) = equatorZ;
        nucdist = [];
        nucelldist = [];
        stack = [];
        for s = 1:length(cellcodes)
%            disp(['s:' num2str(s)]);
            if ~isempty(cellcodes{s})
                curr_nuc = cellcodes{s}.nucdist;
                curr_cell = cellcodes{s}.nucelldist;
                if length(curr_nuc) ~= 360
                    h = -pi:(2*pi/length(curr_nuc)):pi;
                    curr_nuc(end+1) = curr_nuc(1);
                    curr_nuc = interp1(h,curr_nuc,H);
                    curr_nuc(end) = [];
                    curr_cell(end+1) = curr_cell(1);
                    curr_cell = interp1(h,curr_cell,H);
                    curr_cell(end) = [];
                end
                nucdist = [nucdist;curr_nuc];
                nucelldist = [nucelldist;curr_cell];
            end         
        end
        r{i} = nucdist ./ nucelldist;
        
        r{i}(isinf(r{i}(:))) = max(r{i}(~isinf(r{i}(:))));
        
%         if any()
%             1;
%         end
        distratio = [distratio,r{i}(:)];
        cellnucratios = [cellnucratios;cellnucheightratio];
        nucbottomslices = [nucbottomslices;nucbottomslice];
    end
end

rad_ratio = distratio';


end
