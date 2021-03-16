function image = ml_tclread( filename)
% IMAGE = ML_TCLREAD( FILENAME)
% Loads a TCL (.dat) file.

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

% Written by Meel Velliste

f = fopen( filename, 'r', 'ieee-be');
if( f == -1); error(['Could not open file: ' filename]); end;
  
fseek( f, 0, 'bof');
BytesPerPixel = fread( f, 1, 'uint16');
ImageWidth = fread( f, 1, 'uint16');
ImageHeight = fread( f, 1, 'uint16');
SequenceNumber = fread( f, 1, 'uint16');
BitsPerPixel = fread( f, 1, 'uint16');

BEGINNING_OF_DATA = 512;
switch BitsPerPixel
 case 8, PIXEL_TYPE = 'uint8';
 case 16, PIXEL_TYPE = 'uint16';
 otherwise error('BitsPerPixel must be either 8 or 16');
end

fseek( f, BEGINNING_OF_DATA, 'bof');
image = fread( f, [ImageWidth,ImageHeight], PIXEL_TYPE)';

fclose( f);
