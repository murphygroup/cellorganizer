function ml_save3di( img, filename)
% ML_SAVE3DI( IMAGE, FILENAME)
% Saves the matrix IMAGE in the 3D Image file format (.3di)
% which is consists of a simple header (X-dimension, Y-dimension,
% Z-dimension, Bytes per pixel, followed by pixel values).
% Useful for passing an image to a standalone program.

% Copyright (C) 2006  Murphy Lab
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

% Meel Velliste (Spring 2000)

ImageInfo = whos('img');
switch ImageInfo.class
case 'uint8',
    BYTES_PER_PIXEL = 1;
    DATA_TYPE = 'uint8';
    %disp('uint8');
otherwise
    BYTES_PER_PIXEL = 2;
    DATA_TYPE = 'uint16';
end;

f = fopen( filename, 'wb');
if( f == -1) fprintf( 2, 'Cannot open file for saving!\n'); return; end;

s = size( img);
dimensionality = length(s);
if( dimensionality == 2)
    Header = [ s, 1, BYTES_PER_PIXEL];
else
    Header = [ s, BYTES_PER_PIXEL];
end;

fwrite( f, Header, 'ulong');
fwrite( f, img, DATA_TYPE);

fclose( f);
