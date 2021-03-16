function beta = objpos_run( ...
    dnaImagesDirectoryPath, ...
    cellImagesDirectoryPath, ...
    protImagesDirectoryPath, savepath, param )
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

beta = [];
if nargin == 4
    param = [];
    verbose = false;
    debug = false;
elseif nargin > 5
    error('CellOrganizer: Wrong number of input arguments');
else
    try
        verbose = param.verbose;
        if ~islogical( verbose )
            verbose = false;
        end
    catch
        verbose = false;
    end
    
    try
        debug = param.debug;
        if ~islogical( debug )
            debug = false;
        end
    catch
        debug = false;
    end
end

if ~exist(savepath,'dir')
    mkdir(savepath)
end

datapath = [ pwd filesep 'temp/protein_objects_gaussian/object_gaussians'];
%datapath = './inter_results/protein_objects_gaussian/object_gaussians';
gaussobjpath = datapath;

%D. Sullivan 2/22/13
%changed from being hard coded downsample to param based downsampling
% downsample = param.downsampling;
% downsample = [5 5 1];

%icaoberg 15/5/2013
dna_image_files = ml_ls( dnaImagesDirectoryPath );
cell_image_files = ml_ls( cellImagesDirectoryPath );
prot_image_files = ml_ls( protImagesDirectoryPath );
try
    masks_image_files = ml_ls( param.masks );
catch
    mask_image_files = '';
end

if ~exist( [savepath filesep 'beta_all.mat'] )
    X = [];
    Y = [];
    
    imageID = [];
    for i = 1:1:length(cell_image_files)
        try
            disp(['Image ' num2str(i)])
            
            param.cellnum = i;
            %icaoberg 15/5/2013
            try
                dna_image_file = dna_image_files{i};
            catch
                dna_image_file = '';
            end
            
            cell_image_file = cell_image_files{i};
            prot_image_file = prot_image_files{i};
            
            try
                mask_image_file = masks_image_files{i};
            catch
                mask_image_file = [];
            end

            % downsample = [5 5 1];
            %D. Sullivan 2/22/13
            %changed from being hard coded downsample to param based downsampling
            downsample = param.downsampling;
            
            %icaoberg 5/16/2013
            display = false;
            [nucbodyimg, cellbodyimg] = seg_cell_and_nucleus( ...
                dna_image_file, ...
                cell_image_file, ...
                prot_image_file, ...
                mask_image_file, ...
                downsample, display, param );
            param = rmfield( param, 'cellnum' );
            
            nucEdge = bwperim(nucbodyimg,18);
            cellEdge = bwperim(cellbodyimg,18);
            
            % Load object centres
            file = ml_dir( [gaussobjpath filesep '*gaussobjs.mat'] );
            load([gaussobjpath filesep file{1}]);
            objcentres = zeros(size(cellEdge));
            mix = mixes{i};
            offset = offsets{i};
            for p = 1:length(mix)
                %disp(['p:' num2str(p)])
                for t = 1:mix(p).ncentres
                    %disp(['t:' num2str(t)]);
                    ecof = mix(p).centres(t,:) + offset(p,:);
                    ecof = round([ecof(1:2)/5,ecof(3)]);
                    ecof(ecof==0) = 1;
                    ecof(ecof>=200) = 200;
                    objcentres(ecof(1),ecof(2),ecof(3)) = 1;
                end
            end
            objcentres(~cellbodyimg) = 0;
            
            disp('Train position model');
            [distcodes,coords,angles] = ...
                ml_celldistcode(cellEdge,nucEdge,objcentres,1,'all');
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
            save([savepath '_Beta.mat'],'beta')
        catch
            disp(['Ignoring image:' num2str(i)]);
        end
    end
    
    beta = ml_logreg(X,Y)';
    save([savepath filesep 'beta_all.mat'],'beta','X','Y','imageID')
else
    load( [savepath filesep 'beta_all.mat'] );
end
