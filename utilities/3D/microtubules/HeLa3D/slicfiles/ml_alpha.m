function [alphas, signatures, signames] = ml_alpha( input)

% [alphas, signatures, signames] = ml_alpha( input)
%
% Get the alpha signatures of the pointset, where INPUT is either a
% list of 3D coordinates (i.e. a 3xN matrix) or a pathname to a
% pre-computed alpha-shape file.

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

% Run the standalone program that reads the alpha shape files and
% extracts relevant info
getsig = '/home/velliste/matlab/alpha/getsig';
outputfile = ml_fopentemp('/imaging/tmp');
command = [getsig ' ' input ' ' outputfile];
[status, stdout] = unix( command);
if( status ~= 0)
    status, stdout, error('Error in reading the alpha shape files');
end

% Define the types of variable used for signature values
float_type = 'float32';
int_type = 'int32';

% open the alpha output file
f = fopen( outputfile,'r');
if( f == -1) error('Could not open alpha output file for reading'); end

% read the number of alpha values (i.e. length of alpha signature)
[no_of_alphas, count] = fread( f, 1, 'int32');
if( count ~= 1) error('Could not read from alpha output file'); end
no_of_alphas = double( no_of_alphas);

% read the number of float signatures
[no_of_float_sigs, count] = fread( f, 1, 'int32');
if( count ~= 1) error('Could not read from alpha output file'); end
no_of_float_sigs = double( no_of_float_sigs);

% read the number of float signatures
[no_of_int_sigs, count] = fread( f, 1, 'int32');
if( count ~= 1) error('Could not read from alpha output file'); end
no_of_int_sigs = double( no_of_int_sigs);

% read all the float signatures
float_signatures = zeros( no_of_float_sigs, no_of_alphas);
for i = 1 : no_of_float_sigs
    [signature, count] = fread( f, no_of_alphas, float_type);
    if( count ~= no_of_alphas) error('Could not read from alpha output file'); end
    signature = double( signature);
    float_signatures(i,:) = signature';
end

% read all the float signatures
int_signatures = zeros( no_of_int_sigs, no_of_alphas);
for i = 1 : no_of_int_sigs
    [signature, count] = fread( f, no_of_alphas, int_type);
    if( count ~= no_of_alphas) error('Could not read from alpha output file'); end
    signature = double( signature);
    int_signatures(i,:) = signature';
end

% check that end of file reached (otherwise we must have failed to
% read some data)
[last_byte, count] = fread( f, 1, 'uint8');
if( count > 0) error('Excess data in alpha output file'); end

% Close the alpha output file
fclose( f);
unix(['rm -f ' outputfile]);

% Put all the signatures together
signatures = [float_signatures; int_signatures];
% Eliminate first alpha, because the alpha is a negative value
% and the signatures are all 0. Seems irrelevant
signatures(:,1) = [];
% Alpha itself is the first signature
alphas = signatures(1,:);
signatures = signatures(2:end,:);
% Give names to the signatures
signames = {'sig_volume', 'sig_area_f_boundary', 'sig_area_f_regular', ...
	    'sig_area_f_singular', 'sig_length_e_singular', 'sig_w', ...
	    'sig_betti_0', 'sig_betti_1', 'sig_betti_2', 'sig_num_t', ...
	    'sig_num_f_boundary', 'sig_num_f_singular', ...
	    'sig_num_f_regular', 'sig_num_f_interior', ...
	    'sig_num_e_boundary', 'sig_num_e_singular', ...
	    'sig_num_e_regular', 'sig_num_e_interior', ...
	    'sig_num_v_boundary', 'sig_num_v_singular', ...
	    'sig_num_v_regular', 'sig_num_v_interior' }';


