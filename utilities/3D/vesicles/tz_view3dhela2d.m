function h = tz_view3dhela2d(datainfo,classidx,cellidx,param)
%TZ_VIEW3DHELA2D View 3dhela2d data
%   H = TZ_VIEW3DHELA2D(DATAINFO,CLASSIDX,CELLIDX) shows the three channel
%   image of the 3dhela2d data with class index CLASSIDX and cell index
%   CELLIDX. It also returns a handle of the figure. DATAINFO is a
%   structure including essencial information about 
%   
%   H = TZ_VIEW3DHELA2D(DATADIR,CLASSIDX,CELLIDX,PARAM) allows specifying
%   how to view the data. PARAM has the following fields:
%       'option' - option of viewing the data.
%           'rgb' : Show the three channels
%               'rgboption' - see ML_SYNRGBIMG
%           'rgb2' : Show the three channels. But the cell channel is shown as
%               edge.
%           'cellcode' : Show cell code
%               'cellcode' - A structure of the parameters of 
%                   TZ_SHOWCELLCODE
%           'img' : single channel image
%               'channel' - selected channel
%                   'dna' : dna channel
%                   'cell' : total protein channel
%                   'prot' : protein channel
%             
%   See also TZ_SHOWCELLCODE

%   13-Jun-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

% Copyright (C) 2007  Murphy Lab
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

if nargin < 3
    error('3 or 4 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param,struct('option','rgb','rgboption','isc','scale',1));

classDir = [datainfo.dir.root filesep ...
    datainfo.classes{classidx}];
cellSlicePath = [classDir filesep datainfo.dir.cellSlice ...
    filesep datainfo.prefix num2str(cellidx) ...
    datainfo.postfix.cellSlice];
dnaSlicePath = [classDir filesep datainfo.dir.dnaSlice ...
    filesep datainfo.prefix num2str(cellidx) ...
    datainfo.postfix.dnaSlice];
protSlicePath = [classDir filesep datainfo.dir.protSlice ...
    filesep datainfo.prefix num2str(cellidx) ...
    datainfo.postfix.protSlice];

switch param.option
    case 'img'
        switch param.channel
            case 'dna'
                dnaSlice = load(dnaSlicePath);
                imshow(getfield(dnaSlice,datainfo.var.dnaSlice),[]);
            case 'cell'
                cellSlice = load(cellSlicePath);
                imshow(getfield(cellSlice,datainfo.var.cellSlice),[]);
            case 'prot'
                protSlice = load(protSlicePath);
                imshow(getfield(protSlice,datainfo.var.protSlice),[]);   
        end
    case 'rgb'
        cellSlice = load(cellSlicePath);
        dnaSlice = load(dnaSlicePath);
        protSlice = load(protSlicePath);
        h = ml_synrgbimg( ...
            getfield(dnaSlice,datainfo.var.dnaSlice), ...
            getfield(protSlice,datainfo.var.protSlice), ...
            getfield(cellSlice,datainfo.var.cellSlice), ...
            param.rgboption);
        imshow(h,[]);
    case 'rgb2'
        cellEdge = ml_get3dhela2dinfo(datainfo,struct('classidx',classidx, ...
            'cellidx',cellidx,'infotype','data','datatype','celledge'));
        dnaSlice = load(dnaSlicePath);
        protSlice = load(protSlicePath);
        se = strel('disk',4,4);
        cellEdge = imdilate(cellEdge,se);
        h = ml_synrgbimg( ...
            getfield(dnaSlice,datainfo.var.dnaSlice), ...
            getfield(protSlice,datainfo.var.protSlice), ...
            cellEdge, ...
            param.rgboption);
        imshow(h,[]);
    case 'rgb3'
        cellcodePath = [classDir filesep datainfo.dir.cellcode filesep ...
            datainfo.prefix num2str(cellidx) datainfo.postfix.cellcode];
        cellcode = load(cellcodePath);
        cellBoundary = ml_showpts_2d(cellcode.nucellhitpts,'ln',1);
        dnaSlice = load(dnaSlicePath);
        protSlice = load(protSlicePath);
        dnaSlice =  getfield(dnaSlice,datainfo.var.dnaSlice);
        protSlice = getfield(protSlice,datainfo.var.protSlice);
        cellEdge = ml_obj2img(cellBoundary,size(protSlice));
        if param.scale~=1
            cellEdge = imfill(cellEdge,'hole');
            cellEdge = imresize(cellEdge,param.scale,'nearest');
            cellEdge = bwperim(cellEdge);
            dnaSlice = imresize(dnaSlice,param.scale);
            protSlice = imresize(protSlice,param.scale);
        else
            se = strel('disk',2,2);
            cellEdge = imdilate(cellEdge,se);
        end
        
        h = ml_synrgbimg(dnaSlice,protSlice,cellEdge,param.rgboption);
        imshow(h,[]);    
    case 'cellcode'
        if nargout>0
            h = figure;
        end

        cellcodePath = [classDir filesep datainfo.dir.cellcode filesep ...
            datainfo.prefix num2str(cellidx) datainfo.postfix.cellcode];
        cellcode = load(cellcodePath);
        switch param.cellcode.option
            case {'cntr','hitpts','hitxy'}
                tz_showcellcode(cellcode,param.cellcode.option);
            case 'app'
                tz_showcellcode(cellcode,param.cellcode.option, ...
                    param.cellcode.option2,param.cellcode.order);
            case 'img'
                if exist('h','var')
                    close(h);
                end
                img = tz_showcellcode(cellcode,param.cellcode.option, ...
                    param.cellcode.imgsize);
                h = img;
            otherwise
                error('Unrecognized cellcode option');
        end
end

