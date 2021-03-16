function answer = syn2blender( imgfolder, outputfolder, param )
% SYN2BLENDER Exports a generated instance from CellOrganizer to a format
% that can be read by Blender.
%
% List Of Input Arguments     Descriptions
% -----------------------     ------------
% imgfolder                   a folder with a 3d image you wish to obtain the mesh for
% outputfolder                the path where you wish to save the generated files
%
% Parameter structure description
%
% List Of Input Arguments     Descriptions
% -----------------------     ------------
% downsample  (optional) downsample scale. default is 1.
% verbose     (optional) verbose flag that displays progress
% debug       (optional) flag that displays debugging messages. default is false
%
% List Of Outputs     Descriptions
% ---------------     ------------
% answer              true if it saves object file that can be loaded into blender
%                     specified by savefile to disk, false otherwise

% Ivan E. Cao-Berg
%
% Copyright (C) 2012-2019 Murphy Lab
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

% April 21, 2013 D. Sullivan Allowed for the output directory to be created
%                            instead of crashing. Removed warning message.

%icaoberg june 19, 2012
%step 0: check input arguments
answer = false;

if nargin > 3
    error('Wrong number of input arguments');
end

if nargin == 2
    param = struct([]); else
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
    downsample = param.downsample;
    if ~isnumeric( downsample )
        downsample = 1;
    end
catch
    downsample = 1;
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
    %D. Sullivan 4/22/13 allowed for the output directory to be created
    %instead of crashing. Removed warning message
    mkdir(outputfolder)
    %     if debug
    %         warning('Input argument outputfolder does not exist');
    %     end
    %     return
end

%step 1: find all the generated tiff file in imgfolder
files = ml_dir( [ imgfolder filesep '*.tif' ] );

%step 2: make object files

for index=1:1:length(files)
    file = files{index};
    img = tif2img( [imgfolder filesep file] );
    if index == 1
        %R. Arepally 6/5/13 added shiftvector parameter for im2blender
        %function. Used to center objects at origin.
        shiftvector = [];
    end
    savefile = [ outputfolder filesep file(1:end-3) 'obj' ];
    if verbose
        disp( [ 'Exporting ' file ' to ' savefile ] );
    end
    try
        if length(downsample)==1
            if isfield(param,'patchsample')
                if length(param.patchsample) == 1
                    im2blender( img, savefile, downsample,param.patchsample, [], shiftvector );
                elseif length(param.patchsample)==length(files)
                    im2blender( img, savefile, downsample,param.patchsample(index),[], shiftvector);
                else
                    warning('Unable to generate blender files, invalid patchsample vector specified.');
                    return;
                end
            else
                im2blender( img, savefile, downsample, [], shiftvector);
            end
        elseif length(downsample)>1
            if isfield(param,'patchsample')
                if length(param.patchsample) ==1
                    im2blender( img, savefile, downsample(index),param.patchsample);
                elseif length(param.patchsample)==length(files)
                    im2blender( img, savefile, downsample(index),param.patchsample(index));
                else
                    warning('Unable to generate blender files, invalid patchsample vector specified.');
                    return;
                end
            else
                im2blender( img, savefile, downsample(index));
            end
        else
            warning('Unable to generate blender files, invalid downsampling vector specified.');
            return;
        end
    catch
        disp(['Unable to make ' savefile ]);
    end
end

close(gcf);
answer = true;
end%syn2blender
