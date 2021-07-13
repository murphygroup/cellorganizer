function answer = get_murphylab_image_collections( noninteractive )
% GET_MURPHYLAB_IMAGE_COLLECTIONS Helper function that downloads Murphy
% Lab's image collections used by CellOrganizer for model creation and
% demonstrations
%
% The collections are
% * 3D HeLa cells [2.4 GB]
% * 3D movies of T cells expressing LAT (the zip file is 1.2 GB but it
%   expands to 2.6 GB)

% Author: Ivan E. Cao-Berg
%
% Copyright (C) 2015-2017 Murphy Lab
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

[ path, file, extension ] = fileparts( which(mfilename) );
cd( path )
clear path file extension

disp('Checking if image collections are present.');
fileSep = '/';
if ~exist( '.succesfully_downloaded_images' )
    numdownloaded = 0;
    
    %2D/3D HeLa dataset
    tarball = 'cellorganizer_full_image_collection.zip';
    url = 'http://murphylab.web.cmu.edu/data/Hela/3D/multitiff';
    % this just checks for one of the directories and assumes that the rest are there
    if ~exist( [pwd fileSep '3T3'] ) || ~exist( [pwd fileSep 'HeLa'] )
        if nargin == 0
            go=input('Downloading 2D/3D HeLa collection [3.1Gb], enter Y to continue: ', 's' );
        else
            disp('Downloading 2D/3D HeLa collection [3.1Gb]')
            go = 'Y';
        end
        
        if (upper(go)=='Y')
            urlwrite([url fileSep tarball], tarball);
            
            disp('Extracting files');
            unzip( tarball );
            delete( tarball );
            numdownloaded = numdownloaded + 1;
        end
    end
    
    %4D T cell dataset
    tarball = 'LATFull.tgz';
    url = 'http://murphylab.web.cmu.edu/data/TcellModels/';
    if ~exist( [pwd fileSep 'LAT'] )
        if nargin == 0
            go=input('Downloading 4D T cell collection [1.2Gb], enter Y to continue: ','s');
        else
            disp('Downloading 4D T cell collection [1.2Gb]')
            go = 'Y';
        end
        
        if ( upper(go) == 'Y' )
            urlwrite([url fileSep tarball], tarball);
            
            disp('Extracting files');
            !tar -xvf LATFull.tgz
            delete( tarball );
            movefile('./LATFull', './LAT' );
            numdownloaded = numdownloaded + 1;
        end
    end
    
    if numdownloaded == 2
        fclose(fopen('.succesfully_downloaded_images', 'w'));
    end
else
    disp('Image collections already present. Skipping download.');
end
