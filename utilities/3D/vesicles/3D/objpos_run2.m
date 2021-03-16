function beta = objpos_run2( tempfilepath, param )
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
% June 12, 2013 D. Sullivan refactored to run with per-cell precomputed
%                           positions
% July 25, 2013 G. Johnson  Remove NANs from object positions so that
%                           'beta' is always legal
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

%D. Sullivan 6/12/13 don't need this anymore due to refactoring
%D. Sullivan 2/22/13
%changed from being hard coded downsample to param based downsampling
% downsample = param.downsampling;
% downsample = [5 5 1];

%icaoberg 15/5/2013
% dna_image_files = ml_ls( dnaImagesDirectoryPath );
% cell_image_files = ml_ls( cellImagesDirectoryPath );
% prot_image_files = ml_ls( protImagesDirectoryPath );
% try
%     masks_image_files = ml_ls( param.masks );
% catch
%     mask_image_files = '';
% end

if ~exist( [tempfilepath filesep 'beta_all.mat'] )
    Xtot = [];
    Ytot = [];
    files = ml_dir( [tempfilepath filesep '*_Beta.mat'] );
    imageID = [];
    for i = 1:1:length(files)
%         try
            disp(['Image ' num2str(i)])
            load([tempfilepath filesep files{i}]);
            
            %grj 7/23/13
            rminds = any(isnan(X),2) | any(isnan(Y),2);
            
            X(rminds,:) = [];
            Y(rminds,:) = [];
            
            Xtot = [Xtot;X];
            Ytot = [Ytot;Y];
%             betatot(size(betatot,1)+1,:) = ml_logreg(x,distcodes(:,3))';
%             imageID = [imageID;repmat(size(betatot,1)+1,size(x,1),1)];
%         catch
%             disp(['Ignoring image:' num2str(i)]);
%         end
    end
    
    beta = ml_logreg(Xtot,Ytot)';
    save([tempfilepath filesep 'beta_all.mat'],'beta','X','Y','imageID')
else
    load( [tempfilepath filesep 'beta_all.mat'] );
end
