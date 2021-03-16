function answer = syn2projection( imgfolder, outputfolder, param )
% SYN2PROJECTION Makes projections from a set of images synthesized by 
% CellOrganizer
%

% List Of Input Arguments     Descriptions
% -----------------------     ------------
% imgfolder                   a folder of synthesized images by CellOrganizer
% outputfolder                the path where you wish to save the generated files
%
% Parameter structure description
%
% List Of Parameters        Descriptions
% ------------------        ------------
% method                    (optional) either a sum or mean. default is sum
% verbose                   (optional) verbose flag that displays progress
% debug                     (optional) flag that displays debugging messages. default is false
% 
% List Of Outputs     Descriptions
% ---------------     ------------
% answer              true if it saves all projections to disk

% Author: Ivan E. Cao-Berg
%
% Copyright (C) 2012-2016 Murphy Lab
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

%icaoberg june 19, 2012
%step 0: check input arguments
answer = false;

if nargin > 3
    error('Wrong number of input arguments');
end

if nargin == 2
  param = struct([]);
end

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

try
   method = param.method;
catch
   method = 'sum';
end

try
   compression = param.compression;
   if ~strcmpi( compression, 'none' ) && ...
    ~stcmpi( compression, 'lzw' );
      compression = 'none';
   end
catch
   compression = 'none';
end

if isempty( imgfolder )
    if debug
       warning('Input argument imgfolder cannot be empty');
    end
    return
end

if ~exist( imgfolder )
    if debug
       warning('Input argument imgfolder does not exist');
    end
    return
end

if isempty( outputfolder )
    if debug
        warning('Input argument savefile cannot be empty');
    end
    return
end

if ~exist( outputfolder )
    if debug
        warning('Input argument outputfolder does not exist');
    end
    return
end

%step 1: find all the generated tiff file in imgfolder
files = ml_dir( [ imgfolder filesep '*.tif' ] );

%step 2: make object files
for index=1:1:length(files)
   file = files{index};
   image_file = [ imgfolder filesep file ];
   img = tif2img( image_file );
   savefile = [ outputfolder filesep file(1:end-4) '.png' ];
   if verbose
     disp( [ 'Making projection of ' file ] );
   end
   try
       projection = im2projection( img, param );
       imwrite( projection, savefile );
   catch
       disp(['Unable to make ' savefile ]);
   end
end

answer = true;
