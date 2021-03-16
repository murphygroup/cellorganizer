function datainfo = ml_get3dhela2dinfo(data,param)
%ML_GET3DHELA2DINFO get information for 3d hela slice data
%   DATAINFO = ML_GET3DHELA2DINFO(DATA,PARAM) returns information specified
%   by PARAM from DATA which is a structure. PARAM is also a structure and
%   has a field 'infotype', which is a string. Other fields depend on the
%   field 'infotype':
%       'combidx' : indices of the combined data
%           'classidx' - class indices (column vector)
%           'cellidx' - cell indices in the corresponding class (column
%           vector) 
%       'compcell' : if 'compcell' is 1, only complete cells are
%           considered. otherwise only all cells are considered.
%       'classidx' : get class index and cell index
%           'combidx' - a column vector of indices
%       'fileidx' : file indices
%           'classidx' - class indices
%           'cellidx' - cellindices
%       'allcombidx' :  convert indices of complete cells into indices in
%           all cells 
%           'combidx' - 
%       'filepath' : path of a file
%           'classidx' - class index
%           'cellidx' - cell index
%           'datatype' - data types, including 'cellSlice', 'dnaSlice',
%               'protSlice', 'nucedge', 'nucbody', 'celledge', 'cellbody',
%               'cellcode', 'object', 'cellfeat','nuctex','objfeat'
%       'orgfilepath' : orginal file path
%           'classidx - class index
%           'cellidx' - cell index
%           'datatype' - data types,including 'cellSlice', 'dnaSlice',
%               'protSlice', 'mask'
%       'data' : data
%           Other fields are same as those for 'filepath'
%       'rawdata' : original data
%           Other fields are the same as those for 'orgfilepath'
%
%   See also

%   26-Jun-2006 Initial write T. Zhao
%   Copyright (c) 2006 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 2
    error('Exactly 2 arguments are required');
end

param = ml_initparam(param,struct('compcell',0));

if(param.compcell==1)
    for i=1:length(data.incompleteCellIndices)
        cellNumbers(i) = data.cellNumbers(i) - ...
            length(data.incompleteCellIndices{i});
    end
else
    cellNumbers = data.cellNumbers;
end
offsets = [0 cumsum(cellNumbers)];

switch(param.infotype)
    case 'combidx'    
        offsets = offsets';
        datainfo = offsets(param.classidx)+param.cellidx;
    case 'classidx'
        idxmatrix = repmat(param.combidx,1,length(offsets));
        offsetMatrix = repmat(offsets,length(param.combidx),1);
        diffmatrix = idxmatrix - offsetMatrix;
        diffmatrix(diffmatrix<=0) = Inf;
        [datainfo(:,2),datainfo(:,1)] = min(diffmatrix,[],2);
    case 'fileidx'
        if(param.compcell~=1)
            datainfo = param.cellidx;
        else
            for i=1:length(param.classidx)
                incompidx = data.incompleteCellIndices{ ...
                    param.classidx(i)};
                if isempty(incompidx)
                    datainfo(i,1) = param.cellidx(i);
                else
                    datainfo(i,1) = param.cellidx(i)+ ...
                        sum(incompidx<=param.cellidx(i));
                end                
            end            
        end
    case 'filepath'
        for i=1:length(param.classidx)
            datainfo{i} = [data.dir.root filesep ...
                data.classes{param.classidx(i)} ...
                filesep getfield(data.dir,param.datatype) filesep ...
                data.prefix num2str(param.cellidx(i)) ...
                getfield(data.postfix,param.datatype)];
            if isfield(data.var,param.datatype)
                datainfo{i} = [datainfo{i} '.mat'];
            end
        end
    case 'orgfilepath'
        
        for i=1:length(param.classidx)
            datadir = [data.dir.orginal filesep ...
                data.classes{param.classidx(i)} filesep ...
                data.prefix num2str(param.cellidx(i)) ...            
                filesep getfield(data.dir,param.datatype)];
        
            if strcmp(param.datatype,'mask')
                cropfiles = ml_dir([datadir filesep '*.tif*']);
                datainfo{i} = [datadir filesep cropfiles{1}];
            else
                param2 = param;
                param2.infotype = 'filepath';               
            
                filepath = ml_get3dhela2dinfo(data,param2);
                for i=1:length(param.classidx)
                    tmp = load(filepath{i});
                    datainfo{i} = [data.dir.orginal filesep ...
                        data.classes{param.classidx(i)} filesep ...
                        data.prefix num2str(param.cellidx(i)) ...            
                        filesep getfield(data.dir,param.datatype) filesep ...
                        tmp.filename];               
                end
            end
        end
    case 'data'
        switch param.datatype
            case 'nucobj'
                param2 = param;
                param2.datatype = 'nucbody';
                nucbody = ml_get3dhela2dinfo(data,param2);
                param2.datatype = 'dnaSlice';
                nucimg = ml_get3dhela2dinfo(data,param2);
                idx = sub2ind(data.imageSize, ...
                    nucbody(:,1),nucbody(:,2));
                datainfo = [nucbody double(nucimg(idx))];
            otherwise
                param2 = param;
                param2.infotype = 'filepath';
                filepath = ml_get3dhela2dinfo(data,param2);
                if length(filepath)==1
                    loaddata = load(filepath{1});
                    datainfo = getfield( ...
                        loaddata,getfield(data.var,param.datatype));
                else
                    for i=1:length(filepath)
                        loaddata = load(filepath{i});
                        datainfo{i} = getfield( ...
                            loaddata,getfield(data.var,param.datatype));
                    end
                end
        end
    case 'rawdata'
        param2 = param;
        param2.infotype = 'orgfilepath';
        filepath = ml_get3dhela2dinfo(data,param2);
        datainfo = ml_readimage(filepath{1});
    case 'allcombidx'
        cellidx = cell2mat(data.incompleteCellIndices)';
        classidx = [];
        for i=1:length(data.incompleteCellIndices)
            classidx = [classidx;i+zeros( ...
                length(data.incompleteCellIndices{i}),1)];
        end
        incompidx = ml_get3dhela2dinfo(data,struct('infotype','combidx', ...
                'compcell',0,'classidx',classidx,'cellidx',cellidx));
        for i=1:length(param.combidx)
            datainfo(i,1) = ...
                param.combidx(i)+sum(param.combidx(i)>=incompidx);
        end       
    otherwise
        error('Unrecognized datainfomation type.');
end
