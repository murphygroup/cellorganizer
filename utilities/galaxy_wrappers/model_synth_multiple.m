function answer=model_synth_multiple(SYNTHESISFLAG, NUMIMGS, ...
    PREFIX, COMPRESSION, SAMPLINGMETHOD, VERBOSE, ...
    DEBUG, INDEXED, BLENDER, SBML )
% Synthesize instances from multiple trained models saved in multiple formats.

% Xin Lu (xlu2@andrew.cmu.edu)
%
% Copyright (C) 2017 Murphy Lab
% Computational Biology Department
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK
WORKING_DIRECTORY = pwd;
disp( WORKING_DIRECTORY );

files=dir([WORKING_DIRECTORY filesep '*.mat']);

models=sort_nat({files.name});

for ii=1:numel(models)
    models{ii}=[WORKING_DIRECTORY filesep models{ii}];
end

options.targetDirectory = pwd;
options.synthesis = SYNTHESISFLAG;
options.numberOfSynthesizedImages = NUMIMGS;
options.prefix = 'images';
options.image.compression = COMPRESSION;
options.verbose = str2bool(VERBOSE);
options.debug = str2bool(DEBUG);
options.output.tifimages = true;
options.output.indexedimage = str2bool(INDEXED);
options.output.blenderfile = str2bool(BLENDER);
options.output.SBMLSpatial = str2bool(SBML);
options.output.OMETIFF = true;
options.output.blender.downsample = 5;

answer = slml2img( models, options );

system(['find . -type f -name "cell*.ome.tif" -exec mv -v {} . \;']);
system(['find . -type f -name "cell*.xml" -exec mv -v {} . \;']);
system(['find . -type f -name "*.obj" -exec mv -v {} . \;']);

if ~exist( './sbmlspatial' )
	mkdir( './sbmlspatial' );
end

if ~exist( './ometiff' )
	mkdir( './ometiff' );
end

if ~exist( './objs' )
	mkdir( './objs' );
end

system('mv -v *.ome.tif ./ometiff');
system('mv -v *.xml ./sbmlspatial');
system('mv -v *.obj ./objs');

end

function bool=str2bool(str)
if isa( str, 'char' )
    if strcmpi(str,'false')
        bool = false;
    end
    if strcmpi(str,'true')
        bool = true;
    end
else
    bool = str;
end
end