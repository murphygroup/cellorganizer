function answer = slml2report( varargin )
% SLML2REPORT Generate a report comparing two generative models
%
% List Of Input Arguments  Descriptions
% -----------------------  ------------
% model1                   A generative model filename
% model2                   A generative model filename
% options
%
% Example
% > filename1 = '/path/to/model/model1.mat';
% > filename2 = '/path/to/model/model2.mat';
% answer = slml2report( filename1, filename2 );

% Author: Robert F. Murphy
%
% Copyright (C) 2013-2019 Murphy Lab
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
% For additional information visit http://www.cellorganizer.org or
% send email to cellorganizer@compbio.cmu.edu

answer = false;

%parse vargargin
if isdeployed()
       
    filename = is_deployed(varargin{1});
    load(filename);
    
else
    if length(varargin) == 2
        model1_filename = varargin{1};
        model2_filename = varargin{2};
        options = {};
    elseif length(varargin) == 3
        model1_filename = varargin{1};
        model2_filename = varargin{2};
        options = varargin{3};
    else
        warning('Wrong number of input arguments');
        return
    end
end

% init options using ml_initparam
options = ml_initparam(options, struct('verbose', true, 'includenuclear', true, ...
    'includecell', true, 'includeprot', true));

try
    load( model1_filename, 'model' ); model1 = model; clear model;
    load( model2_filename, 'model' ); model2 = model; clear model;
catch err
    warning('Unable to open or read model files')
    getReport(err)
    return;
end

nonuclear = ~isfield(model1,'nuclearShapeModel') | ~isfield(model2,'nuclearShapeModel');
if ~nonuclear
    if ~checkmodelclassandtype(model1.nuclearShapeModel,model2.nuclearShapeModel,'Nuclear')
        return
    end
else
    options.includenuclear = false;
end

nocell = ~isfield(model1,'cellShapeModel') | ~isfield(model2,'cellShapeModel');
if ~nocell
    if ~checkmodelclassandtype(model1.cellShapeModel,model2.cellShapeModel,'Cell')
        return
    end
else
    options.includecell = false;
end

noprotein = ~isfield(model1,'proteinModel') | ~isfield(model2,'proteinModel');
if ~noprotein
    if ~checkmodelclassandtype(model1.proteinModel,model2.proteinModel,'Protein')
        return
    end
else
    options.includeprot = false;
end

models{1} = model1;
models{2} = model2;
classlabels = strcat(model1.name,model2.name);
%if ~noprotein
%    for i=1:2
%      try
%          classlabels = strcat(classlabels, [models{i}.proteinModel.class int2str(i) ';']);
%      catch
%          classlabels = strcat(classlabels, [models{i}.proteinModel.type int2str(i) ';']);
%      end
%      classlabels = classlabels(1:end);
%    end
%end
if exist('./report', 'dir')
    rmdir('./report', 's')
end
mkdir('./report');

fileID = fopen('./index.html', 'w');
html_init(fileID);
header2html(fileID, 'Filenames');
text2html( fileID, ['model1_filename = ' model1_filename] );
fprintf( fileID, ['model2_filename = ' model2_filename] );

models2report_v2(models, options, classlabels, fileID);

html_close(fileID);
fclose(fileID);

if exist([pwd filesep 'index.html'], 'file')
    movefile([pwd filesep 'index.html'], [pwd filesep 'report']);
end

image_info = dir('*.png');
image_names = {image_info.name};
for i = 1:length(image_names)
    movefile(image_names{i}, ['report' filesep image_names{i}]);
end

close all
answer = true;
end%slml2report

function match = checkmodelclassandtype(model1,model2,component)
try
    class_same = strcmp(model1.class, model2.class);
    if ~class_same
        warning([component ' shape model class is not the same.']);
    end
    type_same = strcmp(model1.type, model2.type);
    if ~type_same
        warning([component ' shape model type is not the same.']);
    end
    match = class_same && type_same;
    return
catch
    match = false;
end
end%checkmodelclassandtype