function answer = demo3D30()
% demo3D30
%
% This demo illustrates different ways to sample from points in a
% diffeomorphic model.
%
% Input 
% -----
% * a valid CellOrganizer model file
%
% Output
% ------
% * a random walk

% Devin Sullivan
%
% Copyright (C) 2012-2017 Murphy Lab
% Computational Biology Department
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
if ~isdeployed()
  current_path = which(mfilename);
  [current_path, filename, extension] = fileparts( current_path );
  cd(current_path);
end

disp( 'demo3D30' );
warning('This demo is deprecated. The demo will be removed in future versions of CellOrganizer');
disp( 'The estimated running time is 2 hours 30 minutes. Please wait.' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scriptpath = pwd;
[path,scriptname,ext]=fileparts(mfilename);

load('../../../models/3D/diffeomorphic/helamodel_9_16_15.mat')

p = model.cellShapeModel.positions;

rminds = any(isnan(p),2);

pos = p(~rminds,1:2);

c = 1;
gmmfitfile = 'gmmfit.mat';
if ~exist(gmmfitfile, 'file')
    objgmm = {};
    while c < 5
        tryagain = 0;
        while tryagain < 2
            try
                objgmm{c} = gmdistribution.fit(pos,c, 'replicates', 200);
                break
            catch
                tryagain = tryagain + 1;
            end
        end
        c = c+1;
    end
    save(gmmfitfile, 'objgmm')
else
    load(gmmfitfile)
end

[~, aicind] = min(cellfun(@(x) x.AIC, objgmm));
[~, bicind] = min(cellfun(@(x) x.BIC, objgmm));

gmm_aic = objgmm{aicind};
gmm_bic = objgmm{bicind};

tri = delaunay(pos);

linx = linspace(min(pos(:,1)), max(pos(:,1)), 100);
liny = linspace(min(pos(:,2)), max(pos(:,2)), 100);

[x,y] = meshgrid(linx,liny);

for i = 1:length(x(:))
    gmm_pdf(i) = gmm_aic.pdf([x(i), y(i)]);
end

figure('color', 'w')
axis('tight')

scatter(x(:), y(:), 300, gmm_pdf, '.')
hold on

scatter(pos(:,1), pos(:,2), 300, '.', 'k')
triplot(tri, pos(:,1), pos(:,2), 'k')
saveas(gcf, 'gmmdemo.png', 'png')

positions = nan(size(model.cellShapeModel.positions,1),2);
positions(~rminds,:) = pos;

embed_model = embed_distance_matrix([], struct('positions', model.cellShapeModel.positions(:,1:2)));

model.cellShapeModel.positions = embed_model.positions;
model.nuclearShapeModel.positions = embed_model.positions;
model.cellShapeModel.tessellation = embed_model.tessellation;
model.nuclearShapeModel.tessellation = embed_model.tessellation;

numimgs = 5;
for i = 1:numimgs
    %sample from the gmm ANDed to the convex hull of the triangulation
    simplex_idx = nan;
    while isnan(simplex_idx)
        point = gmm_aic.random;
        simplex_idx = tsearchn(model.cellShapeModel.positions,model.cellShapeModel.tessellation,point);
    end
    points(i,:) = point;
    
    param.position = point;
    
    
    savedir = [scriptpath filesep '_img' num2str(i)];
    if ~exist(savedir, 'dir')
        mkdir(savedir)
    end
    cd(savedir);
    
    param.targetDirectory = [scriptpath filesep scriptname filesep 'img' num2str(i)];
    
    imgs{i} = model2img({model}, param);
end

answer = true;
end
