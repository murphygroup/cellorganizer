function answer = demo3D39()
% demo3D39
%
% This demo illustrates how to sample uniformly at random from a
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

disp('demo3D39');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.targetDirectory = pwd;
scriptpath = pwd;
scriptname = mfilename;

model_filename = '../../../models/3D/diffeomorphic/helamodel_9_16_15.mat';
if exist( model_filename )
    load(model_filename )
    
    cellinds = find(all(~isnan(model.cellShapeModel.positions),2));
    positions = model.cellShapeModel.positions(cellinds,:);
    
    minpos = min(positions,[],1);
    maxpos = max(positions,[],1);
    
    numimgs = 5;
    for i = 1:numimgs %this may take a few seconds, and of course will take long the greater the dimensionality of the space
        %sample from the uniform distribution ANDed to the convex hull of the triangulation
        simplex_idx = nan;
        
        while isnan(simplex_idx)
            point = rand([1, size(minpos,2)]).*(maxpos - minpos) + minpos;
            simplex_idx = tsearchn( model.cellShapeModel.positions, ...
                model.cellShapeModel.tessellation, point );
        end
        
        points(i,:) = point;
        
        options.position = point;
        savedir = [scriptpath filesep '_img' num2str(i)];
        if ~exist(savedir, 'dir')
            mkdir(savedir)
        end
        
        cd(savedir);
        
        options.targetDirectory = [scriptname filesep 'img' num2str(i)];
        options.synthesis.flag = 'framework';
        
        imgs{i} = model2img( {model}, options );
    end
    
    answer = true;
else
    warning(['Model filename ' model_filename ' does not exist.']);
    answer = false;
    return
end