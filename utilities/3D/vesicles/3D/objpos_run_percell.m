function [beta, X, Y, imageID] = objpos_run_percell(segdna, segcell, mixes, offsets, options)
% Learns the protein object position model

% Tao Peng
%
% Copyright (C) 2012-2013 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% April 17, 2012 I. Cao-Berg Added parameter structure to the method and fixed a
%                     bug where the method complained because of a nonexisting
%                     preallocated beta
% Feb 22, 2013 D. Sullivan Changed downsampling to be defined by vector
%                          rather than hard coded
% May 15, 2013 I. Cao-Berg Updated method to support wildcards
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

if ~exist('options', 'var')
    options = [];
end

options = ml_initparam(options, struct('verbose', false, 'debug', false));

beta = [];

options = ml_initparam(options, struct('use_geodesic_distance', false));


%D. Sullivan 6/10/13 gotta assume that this is not necessary, but need to
%figure out which is the extra line

%icaoberg 15/5/2013
% dna_image_files = ml_ls( dnaImagesDirectoryPath );
% cell_image_files = ml_ls( cellImagesDirectoryPath );
% prot_image_files = ml_ls( protImagesDirectoryPath );
% try
%     masks_image_files = ml_ls( param.masks );
% catch
%     mask_image_files = '';
% end

% if ~exist( [savepath filesep num2str(currfile) '_Beta.mat'] )
    X = [];
    Y = [];
    
    imageID = [];
%     for i = 1:1:length(cell_image_files)
%         try
           
            %icaoberg 15/5/2013
%             try
%                 dna_image_file = dna_image_files{i};
%             catch
%                 dna_image_file = '';
%             end
            
%             cell_image_file = cell_image_files{i};
%             prot_image_file = prot_image_files{i};
%             
%             try
%                 mask_image_file = masks_image_files{i};
%             catch
%                 mask_image_file = [];
%             end

            % downsample = [5 5 1];
            %D. Sullivan 2/22/13
            %changed from being hard coded downsample to param based downsampling
            
            %D. Sullivan 6/10/13 removed this call. unnecessary. 
%             [nucbodyimg, cellbodyimg] = seg_cell_and_nucleus( ...
%                 dna_image_file, ...
%                 cell_image_file, ...
%                 prot_image_file, ...
%                 mask_image_file, ...
%                 downsample, display, param );
%             load([param.preprocessingFolder filesep,...
%                 'cell' num2str(currfile)],'segdna','segcell', 'bot_slice'); 
%             
%             segdna = padarray(segdna, [0,0,bot_slice-1], 'pre');
%             segcell = padarray(segcell, [0,0,bot_slice-1], 'pre');
%             
%             clear segcell;
%             clear segdna;
            

            
            nucEdge = bwperim(segdna,18);
            cellEdge = bwperim(segcell,18);
            obj_center_image = zeros(size(cellEdge));
            
            
%             load([param.savefitdir filesep 'gaussobjs_',...
%                 num2str(currfile) '.mat']);
            
            
            mix = mixes{1};
            offset = offsets{1};
            for p = 1:length(mix)
                
                %get the object centers from the gmm
                objcenters = mix(p).centres + repmat(offset(p,:), [mix(p).ncentres, 1]);
                %rescale the object positions to the downsampling
%                 objcenters = objcenters./repmat(options.model.downsampling, [mix(p).ncentres,1]);
                
                objcenters = round(objcenters);
                
                %remove any objects that are outside of the bounds of the image
                objcenters(any(objcenters > repmat(size(cellEdge), [mix(p).ncentres,1])),:) = [];
                objcenters(any(objcenters <= 0,2),:) = [];
                
                %place each object in the image
                for t = 1:size(objcenters)
                    obj_center_image(objcenters(t,1), objcenters(t,2), objcenters(t,3)) = 1;
                end
                
            end
            obj_center_image(~segcell) = 0;
            
            disp('Train position model');
            [distcodes,coords,angles] = ...
                ml_celldistcode(cellEdge,nucEdge,obj_center_image,1,'all', options.use_geodesic_distance);
            nullidx = find(distcodes(:,1)==0 & distcodes(:,2)==0);
            distcodes(nullidx,:)=[];
            %     coords(nullidx,:)=[];
            angles.theta(nullidx,:)=[];
            angles.phi(nullidx,:)=[];
            normdists = distcodes(:,1)./sum(abs(distcodes(:,1:2)),2);
            
            disp('Mapping positions')
            x = ml_mappos([normdists angles.theta angles.phi]);
            X = [X;x];
            Y = [Y;distcodes(:,3)];
            beta(size(beta,1)+1,:) = ml_logreg(x,distcodes(:,3))';
            imageID = [imageID;repmat(size(beta,1)+1,size(x,1),1)];
%             save([savepath filesep num2str(currfile) '_Beta.mat'],...
%                 'beta','X','Y','imageID')
%         catch
%             disp(['Ignoring image:' num2str(i)]);
%         end
%     end
% else 
%     if param.verbose
%         disp(['Found percell betas for cell ',num2str(currfile),'. Skipping.']);
%     end
% end

%D. Sullivan 6/10/13 should be all we need for the total model below. 
% 
%     beta = ml_logreg(X,Y)';
%     save([savepath filesep 'beta_all.mat'],'beta','X','Y','imageID')
% else
%     load( [savepath filesep 'beta_all.mat'] );
% end
