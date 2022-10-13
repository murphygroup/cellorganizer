function answer = slml2info( varargin )
% SLML2INFO Generate a report from information extracted from a genearative model file
%
% List Of Input Arguments  Descriptions
% -----------------------  ------------
% filenames                List of files
% options                  Options structure
%
% Example
% > filenames = {'/path/to/model/file/model.mat'};
% > answer = slml2info( filenames );

% Author: Ivan E. Cao-Berg, Xin Lu, Robert F. Murphy

% Copyright (C) 2020 Murphy Lab
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

%default answer - leave as is
answer = false;

%parse vargargin
if isdeployed()
    is_deployed(varargin)
    
else
    if length(varargin) == 1
        files = varargin{1};
        options = {};
    elseif length(varargin) == 2
        files = varargin{1};
        options = varargin{2};
    else
        warning('Wrong number of input arguments');
        return
    end
end

if ~iscell( files )
    warning('Input argument files must be a cell-array')
    return
end

if length(files) == 1
    filename = files{1};
else
    for i=1:length(files)
        load(files{i});
        if ~is_tcell_model( model )
            warning('Method only accepts one input model of this class/type');
            return
        end
    end
end
%default value is set to report as to not break functionality in
%CellOrganizer for Galaxy and CellOrganizer for Docker. Do not change this
%value
if ~isdeployed() && nargin == 1
    options.output_directory = [pwd filesep 'report'];
end

if ~isfield( options, 'output_directory' )
    options.output_directory = [pwd filesep 'report'];
end

if ~exist( options.output_directory )
    disp(['Creating output directory ' options.output_directory ]);
    mkdir( options.output_directory )
end

outputfile = [ options.output_directory filesep 'index.html' ];

if exist(outputfile, 'file')
    delete(outputfile)
end

fileID = fopen( outputfile, 'w' );
html_init( fileID );

if exist( 'filename', 'var' )
    if exist(filename, 'file')
        load( filename );
        header2html( fileID,'Filename' );
        text2html( fileID, filename );
    else
        disp(['Filename ' filename ' does not exist.']);
        answer = false;
        return
    end
    %when filename exists and files == 1
    model.model_filename=filename;
    model2info_v2( model, fileID, options );
else
    % when files > 1 and filename ! exist
    filename=files{1};
    load( filename ); %this replaces current options
    header2html( fileID,'Filename' );
    text2html( fileID, filename );
    model.model_filename=files;
    model2info_v2( model, fileID, options );
end

html_close( fileID );
fclose( fileID );

if exist([pwd filesep 'html'])
    movefile([pwd filesep 'html' filesep '*'], options.output_directory );
    rmdir([pwd filesep 'html'], 's' );
end
    
answer = true;
end%slml2info
