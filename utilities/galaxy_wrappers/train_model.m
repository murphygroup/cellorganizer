function answer=train_model( DIMENSIONALITY,MASKCHANNEL,DNACHANNEL,CELLCHANNEL,PROTEINCHANNEL, ...
    TRAINFLAG,MODELNAME,MODELID,DOWNSAMPLEX,DOWNSAMPLEY,DOWNSAMPLEZ, ...
    NUCLEUSTYPE,NUCLEUSNAME,NUCLEUSCLASS,NUCLEUSID, ...
    CELLTYPE,CELLNAME,CELLCLASS,CELLID, ...
    PROTEINTYPE,PROTEINNAME,PROTEINCLASS,PROTEINID,PROTEINCYTONUCLEARFLAG, ...
    DOCUMENTATION,VERBOSE,DEBUG )
% function answer = train_model()

% Train a cellular model in CellOrganizer
%
% Xin Lu (xlu2@andrew.cmu.edu)
%
% Copyright (C) 2017 Murphy Lab
% Computational Biology Department
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.o
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

%% parse varargin
options.train.flag = TRAINFLAG; %nuclear, frame, all %% add a default in the shell script so train.flag is always var 3
options.downsampling = [str2num(DOWNSAMPLEX) str2num(DOWNSAMPLEY) str2num(DOWNSAMPLEZ)];
options.nucleus.type = NUCLEUSTYPE; % medial axis, cylindrical surface or diffeo
options.nucleus.class = NUCLEUSCLASS; % nuc or framework
if strcmp(options.train.flag,'framework')
    options.cell.type = CELLTYPE;
    options.cell.class = CELLCLASS;
end
if strcmp(options.train.flag,'all')
    options.protein.type = PROTEINTYPE;
    options.protein.class = PROTEINCLASS;
    options.protein.cytonuclearflag = PROTEINCYTONUCLEARFLAG;
end
options.verbose = str2bool(VERBOSE);
options.debug = str2bool(DEBUG);
options.model.name = MODELNAME; %same as above, optional
options.model.filename = 'model.mat';
options.model.id = MODELID; %random string%
options.nucleus.name = NUCLEUSNAME; %defualt to empty
options.nucleus.id = NUCLEUSID;
options.cell.name = CELLNAME;
options.cell.id = CELLID;
options.protein.name = PROTEINNAME;
options.protein.id = PROTEINID;
options.documentation = DOCUMENTATION;

if strcmpi(options.model.id, '')
    options.model.id = num2str(now);
end
if strcmpi(options.cell.id, '')
    options.cell.id = num2str(now);
end
if strcmpi(options.nucleus.id, '')
    options.nucleus.id = num2str(now);
end
if strcmpi(options.protein.id, '')
    options.protein.id = num2str(now);
end

%% pairs image data
directory=pwd;
ometiff_path=[directory filesep '*.ome.tif'];
files = dir(ometiff_path);
files = sort_nat({files.name});

cell={};
protein={};
if has_rois(files{1})
    dna=get_list_of_function_handles_from_wildcards(ometiff_path,1);
    if strcmp(options.train.flag,'framework')
        cell=get_list_of_function_handles_from_wildcards(ometiff_path,2);
    end
    if strcmp(options.train.flag,'all')
        protein=get_list_of_function_handles_from_wildcards(ometiff_path,3);
    end
else
    for i = 1:length(files)
        dna{i} = eval(['@() uint8(flipdim(OME_loadchannel([directory filesep files{i}],' num2str(DNACHANNEL) '),3))']);
        if strcmp(options.train.flag,'framework')
            cell{i} = eval(['@() uint8(flipdim(OME_loadchannel([directory filesep files{i}],' num2str(CELLCHANNEL) '),3))']);
        end
        if strcmp(options.train.flag,'all')
            protein{i} = eval(['@() uint8(flipdim(OME_loadchannel([directory filesep files{i}],' num2str(PROTEINCHANNEL) '),3))']);
        end
        masks{i} = eval(['@() uint8(flipdim(OME_loadchannel([directory filesep files{i}],' num2str(MASKCHANNEL) '),3))']);
    end
    options.masks = masks;
end

% get image resolution and dimensionality
% options.model.resolution = [0.049, 0.049, 0.2000];
options.model.resolution = OME_getResolution([directory filesep files{1}])
if ~isempty(options.model.resolution)
    dimensionality = DIMENSIONALITY;
    save( 'options.mat', 'dimensionality', 'options', 'dna', 'cell', 'protein' );
    answer = img2slml(dimensionality, dna, cell, protein, options);
else
    disp('No resolution specified, exiting.');
end
end

function answer=str2bool(str)
if isa( str, 'char' )
    if strcmpi(str,'false')
        answer = false;
    end
    if strcmpi(str,'true')
        answer = true;
    end
else
    bool = str;
end
end