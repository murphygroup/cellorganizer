function image = ml_readimage( filename, tmpdir )
% IMAGE = ML_READIMAGE( FILENAME,TMPDIR)
%
% Reads a 2D image in any format, including TCL files.
% If FILENAME is an emtpy string, then an empty matrix
% is returned.

% Copyright (C) 2006-2013 Murphy Lab
% Carnegie Mellon University
%
% May 1, 2013 I. Cao-Berg Updated method so that it will try to open the
% image using Bio-Formats first, and if errors occurs it will try to open
% the file with imread
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

% ml_readimage written by Meel Velliste
%G. Johnson 10/23/13 Added support for function handles
%T. Majarian 7/26/2016 Added support for OME Tif + channel label, input can
%be in cell form {filepath, channel_name}

if(nargin<2)
    tmpdir = tempdir;
end

if( tmpdir(end) ~= filesep)
    tmpdir(end+1) = filesep;
end

if strcmpi(class(filename), 'function_handle')
    image = feval(filename);
    return;
elseif( isempty( filename) || ~isa( filename, 'char' ) )
    image = [];
    return;
%%%% support for OME tif
elseif strfind(filename, ',')
    splitpath = strsplit(filename, ',');
    image_path = splitpath{1};
    if exist(image_path,'file')
        image = OME_loadchannel(image_path,splitpath{2});
        return;
    end
end

L = length( filename);
extension = filename(L-2:L);

if isempty( extension )
    image = [];
else
    switch( extension )
        case '.gz'
            [f,tempfilename] = ml_fopentemp(tmpdir);
            tempfullpath = [tmpdir tempfilename];
            fclose(f);
            unix(['gunzip -c ' filename ' > ' tempfullpath]);
            image = imread( tempfullpath);
            unix(['rm ' tempfullpath]);
        case 'bz2'
            [f,tempfilename] = ml_fopentemp(tmpdir);
            tempfullpath = [tmpdir tempfilename];
            fclose(f);
            unix(['bunzip2 -c ' filename ' > ' tempfullpath]);
            image = imread( tempfullpath);
            unix(['rm ' tempfullpath]);
        case 'dat'
            image = ml_tclread( filename);
        otherwise
            try
                %it checks if the file is a multitiff, if so it tries to
                %open the file using tif2img rather than imread
                info = imfinfo( filename );
                if length(info) > 1 && strcmpi( info(1).Format, 'tif' )
                    image = tif2img( filename );
                else
                    image = imread( filename );
                end
            catch
                bfimg = bfopen( filename );
                if ~isempty( bfimg )
                    number_of_images = length( bfimg{1} );
                    image = [];
                    
                    for i=1:1:number_of_images
                        if isempty( image )
                            image(:,:,size(image,3) ) = bfimg{1}{i};
                        else
                            image(:,:,size(image,3)+1 ) = bfimg{1}{i};
                        end
                    end
                else
                    image = [];
                end
            end
    end
end