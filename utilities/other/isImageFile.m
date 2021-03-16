function answer = isImageFile( file )
%ISIMAGEFILE True iff file is supported by the CellOrganizer, i.e. tiff, jpg and dat files.

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
%
% Copyright (C) 2008-2012 Murphy Lab
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

answer = false;
if( isTiff( file ) || isDat( file ) )
    answer = true;
end
end%isImageFile

function answer = isTiff( file )
%ISTIFF True iff file exists and is a tiff file, false otherwise.
if( nargin == 1 )
    if( ~exist( file, 'file' ) )
        error( ['CellOrganizer: The file ' file ' does not exist.'] );
    else
        try
            info = imfinfo( file );
            if( strcmpi( info.Format, 'tif' ) )
                answer = true;
            else
                answer = false;
            end
        catch
            if( strcmpi( getFileExtension( file ), 'tif' ) )
                answer = true;
            else
                answer = false;
            end
        end
    end
else
    error( 'CellOrganizer: Wrong number of input arguments.' );
end
end%isTiff

function answer = isPng( file )
%ISPNG True iff file exists and is a png file, false otherwise.

if( nargin == 1 )
    if( ~exist( file, 'file' ) )
        error( ['CellOrganizer: The file ' file ' does not exist.'] );
    else
        try
            info = imfinfo( file );
            if( strcmpi( info.Format, 'jpg' ) || ...
                    strcmpi( info.Format, 'png' ) )
                answer = true;
            else
                answer = false;
            end
        catch
            if( strcmpi( getFileExtension( file ), 'png' ) || ...
                    strcmpi( getFileExtension( file ), 'png' ) )
                answer = true;
            else
                answer = false;
            end
        end
    end
else
    error( 'CellOrganizer: Wrong number of input arguments.' );
end
end%isPng

function answer = isDat( file )
%ISDAT True iff file exists and is a dat file, false otherwise.

if( nargin == 1 )
    if( ~exist( file, 'file' ) )
        error( ['CellOrganizer: The file ' file ' does not exist.'] );
    else
        try
            info = imfinfo( file );
            if( strcmpi( info.Format, 'dat' ) )
                answer = true;
            else
                answer = false;
            end
        catch
            if( strcmpi( getFileExtension( file ), 'dat' ) )
                answer = true;
            else
                answer = false;
            end
        end
    end
else
    error( 'CellOrganizer: Wrong number of input arguments.' );
end
end%isDat
