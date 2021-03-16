function images = getImages( directory, extension )
%GETIMAGES Returns a cell array containing all the images files in the
%directory.

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
% Created: September 17, 2007
% Last Update: March 2, 2008
%
% Copyright (C) 2008 Center for Bioimage Informatics/Murphy Lab
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

images = {};
%getImages( pwd )
if( nargin == 0 )
    files = getFiles( pwd );
    for i = 1:1:length( files )
        if( isImageFile( files{i} ) )
            images{length(images)+1} = files{i};
        end
    end
    % getImages( directory )
elseif( nargin == 1 )
    if( ~isdir( directory ) )
        warning( 'CellOrganizer: Directory nonexistent or invalid' );
        return
    else
        files = getFiles( directory );
        for i = 1:1:length( files )
            if( isImageFile( files{i} ) )
                images{length(images)+1} = files{i};
            end
        end
    end
    % getImages( directory, extension )
elseif( nargin == 2 )
    if( ~isdir( directory ) )
        warning( 'CellOrganizer: Directory nonexistent or invalid' );
        return
    else
        files = getFiles( directory );

        switch lower(extension)
            case {'jpg', 'jpeg'}
                for i = 1:1:length( files )
                    if( isJpg( [directory filesep files{i}] ) )
                        images{length(images)+1} = files{i};
                    end
                end
            case 'bmp'
                for i = 1:1:length( files )
                    if( isBmp( [directory filesep files{i}] ) )
                        images{length(images)+1} = files{i};
                    end
                end
            case 'dat'
                for i = 1:1:length( files )
                    if( isDat( [directory filesep files{i}] ) )
                        images{length(images)+1} = files{i};
                    end
                end
            case {'tif', 'tiff'}
                for i = 1:1:length( files )
                    if( isTiff( [directory filesep files{i}] ) )
                        images{length(images)+1} = files{i};
                    end
                end
            case 'png'
                for i = 1:1:length( files )
                    if( isPng( [directory filesep files{i}] ) )
                        images{length(images)+1} = files{i};
                    end
                end
            case 'gif'
                for i = 1:1:length( files )
                    if( isGif( [directory filesep files{i}] ) )
                        images{length(images)+1} = files{i};
                    end
                end
            case 'raw'
                for i = 1:1:length( files )
                    if( isRaw( [directory filesep files{i}] ) )
                        images{length(images)+1} = files{i};
                    end
                end
            otherwise
                error( 'CellOrganizer: Unknown or unsupported file extension.' );
        end
    end
else
    error( 'CellOrganizer: Wrong number of input arguments' );
end
end%getImages