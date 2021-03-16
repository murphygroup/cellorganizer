function answer = model2slml( model, filename, param )
%MODEL2SLML Reads a trained generative model of subcellular location from
%CellOrganizer and exports it as a valid SLML instance.
%
%Inputs
%model: a matlab structure holding the model
%filename: the output filename
%
%Output
%answer: true if the model was saved to disk successfully, false otherwise
%
%See also SLML2MODEL

% Author: Yue Yu
%
% Copyright (C) 2012 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% 7/11/2012 Y. Yu part of the code is taken from Michelle's code:param2xml.m
% 7/15/2012 Y. Yu add a type attribute to vector node
% 7/16/2012 Y. Yu Fixed a bug when reading 2D model, now it can read models
%              which contain cell array
% 1/18/2013 I. Cao-Berg Fixed names of classes: referenced class names 
%              rather than instance names
% 1/20/2012 I. Cao-Berg Added try/catch statements at making the document
% object model and parsing the xml files
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published1
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

%step 1: check input arguments
if nargin == 2
    param = struct([]);
elseif nargin == 3
    if ~isstruct( param )
        warning( 'Input argument param must be a structure.' );
        param = struct([]);
    end
else
    error( 'Wrong number of input arguments.' );
end

answer = false;
if isempty( model )
    warning( 'Input argument model cannot be empty. Exiting method.' );
    return
end

if ~isstruct( model )
    warning( 'Input argument model has to be a structure. Exiting method' );
    return
end

if isempty( filename )
    warning( 'Input argument filename cannot be empty. Exiting method.' );
    return
end

if ~isa( filename, 'char' )
    warning( 'Input argument filename must a string. Exiting method.' );
    return
end

try
    verbose = param.verbose;
catch
    verbose = true;
end

try
    debug = param.debug;
catch
    debug = false;
end

%if the input is a single model, make it the required format
if isstruct( model )
    try
        modelname = model.name;
    catch
        modelname = 'unknown model name';
    end
    
    lom = {modelname,model};
    cell.lom = lom;
    
    LOC = {'cell1',cell};
    
    model = LOC;
end

%step 2: check the structure holds a valid slml model
if verbose
    disp( 'Constructing DOM model' );
    tic;
end

%the meat and bones of this method. this submethod creates the document
%object model (DOM) that will be saved as an XML file
try
    root_node= read_LOC( model );
    
    if verbose
        toc
    end
    
    return
catch
    warning('Unable to create document object model');
    if verbose
        toc
    end
    
    return
end

%step 3: save model to disk
%even though this step is not trivial, it is very straight forward.
%once the document object model (DOM) is built correctly, this method
%takes care of parsin the object to an XML file. failure is not
%expected at this point with the exception of memory leaks which are
%unlikely but plausible
if verbose
    disp( 'Parsing DOM model to XML file and saving to disk' );
    tic
end

try
    xmlwrite(filename,root_node);
    if verbose
        toc
    end
catch
    warning('Unable to save document object model to XML file');
    
    if verbose
        toc
    end
    
    return
end

%step 4: response to user
%if failure hasn't occurred at this point it is safe to assume that
%everything went according to plan, return a positive message
answer = true;
end

function root_node = read_LOC( LOC, debug )
%READ_LOC Helper method that reads a list of cells into a DOM
%@param LOC the list of cells, cell array {name,cell,name,cell...} each 
%           cell is defined in read_cell function
%@return the root node

%step1 : check input parameters
if nargin < 2
    debug = false; %set default value
end

root_node = [];
if isempty( LOC ) % to check if the input dictionary is empty
    if debug
        error('Error: Input cell array is empty')
    end
end

N = length(LOC)/2;

%initialize the node
[root_node, cell_node] = read_cell(LOC{2});
name_attrnode = root_node.createAttribute('name');
name_attrnode.setNodeValue(LOC{1});
cell_node.setAttributeNode(name_attrnode);
LOC_node = root_node.createElement('listOfCells');
LOC_node.appendChild(cell_node);

for i = 1 : N-1
    [root_node, cell_node] = read_cell(LOC{i*2+2},root_node);
    name_attrnode = root_node.createAttribute('name');
    name_attrnode.setNodeValue(LOC{i*2+1});
    cell_node.setAttributeNode(name_attrnode);
    LOC_node.appendChild(cell_node);
end
root_node.getDocumentElement.appendChild(LOC_node);
root_attr = root_node.getDocumentElement;
root_attr.setAttribute('level','1');
root_attr.setAttribute('version','2.0');
end

function [root_node, cell_node] = read_cell(cell, root_node, debug)
% read a cell into a DOM object
% @param cell, a matlab struct, required field:
% lom, the list of models, a cell array,{name,model,name,model...} each element is a model which
% is defined in read_model function
% @return a cell DOM object
if nargin < 2
    debug = false; %set default value
    root_node = [];
end
if nargin == 2
    debug = false; %set default value
end
%step1 : check input parameters
if isempty(cell) % to check if the input dictionary is empty
    if debug
        error('Error: Input dictionary is empty')
    end
end
if ~isstruct(cell) % to check if the input is dictionary
    if debug
        error('Error: Input argument must be a struct')
    end
end
if ~isfield(cell,'lom')
    if debug
        error('Error: a lom field is required')
    end
end

lom = cell.lom;
N = length(lom)/2;
cell = rmfield(cell,'lom');

%icaoberg ignore documentation structure
if isfield( lom{2}, 'documentation' )
    lom{2} = rmfield( lom{2}, 'documentation' );
end

%initialize the node
if isempty(root_node)
    [root_node, model_node] = read_model(lom{2});
else
    [root_node, model_node] = read_model(lom{2},root_node);
end
lom_node = root_node.createElement('listOfModels');
name_attrnode = root_node.createAttribute('name');
name_attrnode.setNodeValue(lom{1});
model_node.setAttributeNode(name_attrnode);
lom_node.appendChild(model_node);
%read all models
for i = 1 : N-1
    [root_node, model_node] = read_model(lom{i*2+2},root_node);
    name_attrnode = root_node.createAttribute('name');
    name_attrnode.setNodeValue(lom{i*2+1});
    model_node.setAttributeNode(name_attrnode);
    lom_node.appendChild(model_node);
end
%set all attribute
cell_node = root_node.createElement('cell');
cell_node.appendChild(lom_node);
names = fieldnames(cell);
NN = length(names);
for i = 1 : NN
    tmp_attrnode = root_node.createAttribute(names{i});
    tmp_attrnode.setNodeValue(getfield(cell,names{i}));
    cell_node.setAttrbuteNode(tmp_attrnode);
end
end

function [root_node, model_node] = read_model(model, root_node,debug)
%  read a model into a DOM object
%  @param model, a matlab struct, all of its field which are also
%  struct will be considered as compartment. all other field will be
%  consider as attribute.
%  @return a model DOM node
if nargin < 2
    debug = false; %set default value
    root_node = [];
end
if nargin == 2
    debug = false; %set default value
end
%step1 : check input parameters
if isempty(model) % to check if the input dictionary is empty
    if debug
        error('Error: Input dictionary is empty')
    end
end
if ~isstruct(model) % to check if the input is dictionary
    if debug
        error('Error: Input argument must be a struct')
    end
end
names = fieldnames(model); % get all field names
N = length(names); % number of fields
loc = {};

for i = 1 : N
    tmp = getfield(model,names{i});
    if isstruct(tmp)
        model = rmfield(model,names{i});
        loc{end+1} = names{i};
        loc{end+1} = tmp;
    end
end

if isempty(root_node)
    [root_node,loc_node] = read_loc(loc);
else
    [root_node,loc_node] = read_loc(loc,root_node);
end

model_node = root_node.createElement('model');
model_node.appendChild(loc_node);

%set the attributes
attrnames = fieldnames(model);
NN = length(attrnames);

for i = 1 : NN
    tmp_attrnode = root_node.createAttribute(attrnames{i});
    tmp_attrnode.setNodeValue(getfield(model,attrnames{i}));
    model_node.setAttributeNode(tmp_attrnode);
end

end

function [root_node,loc_node] = read_loc(loc, root_node,debug)
% read a list of compartment into a DOM object
% @param the list of compartment(loc), cell array, {name,comp,name,comp....}
% each element is a compartment which is defined in read_comp function
% @return a DOM node.

if nargin < 2
    debug = false; %set default value
    root_node = [];
end
if nargin == 2
    debug = false; %set default value
end

%step1 : check input parameters
if isempty(loc) % to check if the input dictionary is empty
    if debug
        error('Error: Input dictionary is empty')
    end
end


N = length(loc)/2; % the number of compartment
%initialize the node
if isempty(root_node)
    [root_node,compart_node] = read_comp(loc{2});
else
    [root_node,compart_node] = read_comp(loc{2},root_node);
end
loc_node = root_node.createElement('listOfCompartment');
name_attr = root_node.createAttribute('name');
name_attr.setNodeValue(loc{1});
compart_node.setAttributeNode(name_attr);
loc_node.appendChild(compart_node);

for i = 1 : N-1
    [root_node,compart_node] = read_comp(loc{i*2+2},root_node);
    name_attr = root_node.createAttribute('name');
    name_attr.setNodeValue(loc{i*2+1});
    compart_node.setAttributeNode(name_attr);
    loc_node.appendChild(compart_node);
end
end

function [root_node,compart_node] = read_comp(compartment, root_node, debug)
% read a compartment into a DOM object
% @param the compartment(compartment), matlab struct, which has four possible field:
% shape, texture,position and frequency.
% each of them is a list of parameters and the format is same as the input
% of read_lop function.if none of those four field exist in the compartment,
% then the whole compartment is considered as one list of parameters and the name is 'shape'.
% @return a DOM node.
if nargin < 2
    debug = false; %set default value
    root_node = [];
end
if nargin == 2
    debug = false; %set default value
end
%step1 : check input parameters
if isempty(compartment) % to check if the input dictionary is empty
    if debug
        error('Error: Input dictionary is empty')
    end
end
if ~isstruct(compartment) % to check if the input is dictionary
    if debug
        error('Error: Input argument must be a dictionary')
    end
end
% if the compartment doesn't have any of those four field, then the whole
% compartment is considered as 'shape'
if ~isfield(compartment,'size') && ~isfield(compartment,'texture') && ...
        ~isfield(compartment,'position') && ~isfield(compartment,'frequency') &&...
        ~isfield(compartment,'positionModel') && ~isfield(compartment,'objectModel')
    if isempty(root_node)
        [root_node,lop_node] = read_lop(compartment,'shape');
    else
        [root_node,lop_node] = read_lop(compartment,'shape',root_node);
    end
    shape_node = root_node.createElement('shape');
    shape_node.appendChild(lop_node);
    
    compart_node = root_node.createElement('compartment');
    compart_node.appendChild(shape_node);
    
    %root_node.getDocumentElement.appendChild(compart_node);
    %xmlwrite('compart_try.xml',root_node);
    return
end
component = {}; % used for store each component,{name, value, name, value...}
%for 3D model
if isfield(compartment,'size') && isstruct(compartment.size)
    component{end + 1} = 'shape';
    component{end + 1} = compartment.size;
    compartment = rmfield(compartment,'size');
end
%for 2D model
if isfield(compartment,'objectModel') && isstruct(compartment.objectModel)
    component{end + 1} = 'shape';
    component{end + 1} = compartment.objectModel;
    compartment = rmfield(compartment,'objectModel');
end
if isfield(compartment,'texture') && isstruct(compartment.texture)
    component{end + 1} = 'texture';
    component{end + 1} = compartment.texture;
    compartment = rmfield(compartment,'texture');
end
%for 3D model
if isfield(compartment,'position') && isstruct(compartment.position)
    component{end + 1} = 'position';
    component{end + 1} = compartment.position;
    compartment = rmfield(compartment,'position');
end
%for 2D model
if isfield(compartment,'positionModel') && isstruct(compartment.positionModel)
    component{end + 1} = 'position';
    component{end + 1} = compartment.positionModel;
    compartment = rmfield(compartment,'positionModel');
end
if isfield(compartment,'frequency') && isstruct(compartment.frequency)
    component{end + 1} = 'frequency';
    component{end + 1} = compartment.frequency;
    compartment = rmfield(compartment,'frequency');
end
N = length(component)/2; % get the number of component

%initialize the node
init_name = component{1};
init_value = component{2};
if isempty(root_node)
    [root_node,init_node] = read_lop(init_value,init_name);
else
    [root_node,init_node] = read_lop(init_value,init_name,root_node);
end
compo_node = root_node.createElement(init_name);
compo_node.appendChild(init_node);

compart_node = root_node.createElement('compartment');
compart_node.appendChild(compo_node);

for i = 1 : N-1
    tmp_name = component{i*2 + 1};
    tmp_value = component{i*2 + 2};
    [root_node,tmp_node] = read_lop(tmp_value,tmp_name,root_node);
    
    compo_node = root_node.createElement(tmp_name);
    compo_node.appendChild(tmp_node);
    compart_node.appendChild(compo_node);
end
% set all other field as attribute
names = fieldnames(compartment); % get all fields names
numofattr = length(names);% the number of attrs
for i = 1 : numofattr
    tmp_attr = getfield(compartment,names{i});
    tmp_attrnode = root_node.createAttribute(names{i});
    tmp_attrnode.setNodeValue(tmp_attr);
    compart_node.setAttributeNode(tmp_attrnode);
end

%root_node.getDocumentElement.appendChild(compart_node);
%xmlwrite('compart_try.xml',root_node);
end

function flat_struct = flatstruct(structure,flat_struct,currname)
% flat the struct or cell, set all field name as structure.attr.attr....  form
% Input: structure, a matlab struct(not necessary during the recursion process)
%        flat_struct, the cell array used for storing information, initial
%        should be null cell array
%        currname, string, the variable name of current variable
% Output: a matlab cell array, {name,value,name,value...}
if iscell(structure) % make it structure
    numoffield = length(structure);
    new_structure = struct;
    for i = 1 : numoffield
        tmp_name = strcat('cellstruct',num2str(i));
        new_structure.(tmp_name) = structure{i};
    end
    flat_struct = flatstruct(new_structure,flat_struct,currname);
end
if ~isstruct(structure) && ~iscell(structure) % the current data is not a struct or cell
    flat_struct(end+1) = {currname};
    flat_struct(end+1) = {structure};
    return
end
if isstruct(structure)
    %the current data is struct
    F = fieldnames(structure);%get all field names
    N = length(F); % get the number of fields
    for i = 1 : N
        tmp_field_value = getfield(structure,F{i});
        flat_struct = flatstruct(tmp_field_value,flat_struct,strcat(currname,strcat('.',F{i})));
    end
end
end

% function flat_struct = flatstruct(structure,flat_struct,currname)
% % flat the struct, set all field name as structure.attr.attr....  form
% % Input: structure, a matlab struct(not necessary during the recursion process)
% %        flat_struct, the cell array used for storing information, initial
% %        should be null cell array
% %        currname, string, the variable name of current variable
% % Output: a matlab cell array, {name,value,name,value...}
%
% if ~isstruct(structure) % the current data is not a struct or cell
%     flat_struct(end+1) = {currname};
%     flat_struct(end+1) = {structure};
%     return
% end
% %the current data is struct
% F = fieldnames(structure);%get all field names
% N = length(F); % get the number of fields
% for i = 1 : N
%     tmp_field_value = getfield(structure,F{i});
%     flat_struct = flatstruct(tmp_field_value,flat_struct,strcat(currname,strcat('.',F{i})));
% end

function [root_node,lop_node] = read_lop(lop,name,root_node,debug)
%READ_LOP Helper method that converts a bunch of parameters into a DOM node.
%     @param name, the name of current lop
%     @param root_node, the root node used for generating other lops,
%     default is false
%     @param a list of parameters(lop), a matlab struct, field-values:
%     {name:string(default none),id:string(default none), comment:string(default none),
%     all other fields are
%     considered as parameter values using same parameter name
%     }
%     for each parameter:
%     could be any structure
%     {name:string(default none),id:string(default none),comment:string(default none)
%     type:string(int, string, double),isNumeric: true/false.
%     class: string(scalar,matrix,vector)
%     all other fields are
%     considered as parameter values }
%     @output an DOM object, the root node

%step1 : check input parameters

%debug flag set to default value: false
if nargin <= 2
    debug = false;
    root_node = [];
end

if nargin <= 3
    debug = false; %set default value
end

if isempty( lop ) % to check if the input dictionary is empty
    if debug
        error('Error: Input dictionary is empty')
    end
end

if ~isstruct(lop) % to check if the input is dictionary
    if debug
        error('Error: Input argument must be a dictionary')
    end
end

if isempty(root_node) % if the root_node is not provided, initialize it
    root_node = com.mathworks.xml.XMLUtils.createDocument('slml');
end

lop_node = root_node.createElement('listOfParameters');
% check if name, id or comment exist
if isfield( lop,'name' )
    lop_name = lop.name;
    lop = rmfield( lop,'name' );
    tmp_node = root_node.createAttribute( 'name' );
    tmp_node.setNodeValue( lop_name );
    lop_node.setAttributeNode( tmp_node );
elseif isfield(lop,'id')
    lop_id = lop.id;
    lop = rmfield(lop,'id');
    tmp_node = root_node.createAttribute('id');
    tmp_node.setNodeValue(lop_id);
    lop_node.setAttributeNode(tmp_node);
elseif isfield(lop,'comment')
    lop_comment = lop.comment;
    lop = rmfield(lop,'comment');
    tmp_node = root_node.createAttribute('comment');
    tmp_node.setNodeValue(lop_comment);
    lop_node.setAttributeNode(tmp_node);
end

%retrieve all parameters, store it as cell array
flat_parameter = flatstruct(lop,{},name);
entry_num = length(flat_parameter);
if mod(entry_num,2) ~= 0
    if debug
        disp('wrong flat structure')
    end
    return
end

numofpairs = entry_num/2;
for i = 0 : numofpairs-1
    curr_name = flat_parameter{i*2+1};
    curr_value = flat_parameter{i*2+2};
    [a,b] = size(curr_value); % get the size of current value
    param_node = root_node.createElement('parameter');
    
    % set all attributes
    name_attr = root_node.createAttribute('name');
    name_attr.setNodeValue(curr_name);
    param_node.setAttributeNode(name_attr);
    % ignore the id and comment for now
    if ischar(curr_value) % the current value is string
        type_attr = root_node.createAttribute('type');
        type_attr.setNodeValue('string');
        param_node.setAttributeNode(type_attr);
        
        isNum_attr = root_node.createAttribute('isNumeric');
        isNum_attr.setNodeValue('false');
        param_node.setAttributeNode(isNum_attr);
        
        class_attr = root_node.createAttribute('class');
        class_attr.setNodeValue('scalar');
        param_node.setAttributeNode(class_attr);
    else
        if floor(curr_value) == curr_value
            type_attr = root_node.createAttribute('type');
            type_attr.setNodeValue('integer');
            param_node.setAttributeNode(type_attr);
        else
            type_attr = root_node.createAttribute('type');
            type_attr.setNodeValue('double');
            param_node.setAttributeNode(type_attr);
        end
        
        isNum_attr = root_node.createAttribute('isNumeric');
        isNum_attr.setNodeValue('true');
        param_node.setAttributeNode(isNum_attr);
        
        if a == 1 && a == b
            class_attr = root_node.createAttribute('class');
            class_attr.setNodeValue('scalar');
            param_node.setAttributeNode(class_attr);
        elseif a > 1 && b > 1
            class_attr = root_node.createAttribute('class');
            class_attr.setNodeValue('matrix');
            param_node.setAttributeNode(class_attr);
        else
            class_attr = root_node.createAttribute('class');
            class_attr.setNodeValue('vector');
            param_node.setAttributeNode(class_attr);
        end
    end
    
    % store the value
    if strcmp(class_attr.getNodeValue,'scalar') % store it as scalar
        scalar_node = root_node.createElement('scalar');
        text_node = root_node.createTextNode(num2str(curr_value));
        scalar_node.appendChild(text_node);
        param_node.appendChild(scalar_node);
    elseif strcmp(class_attr.getNodeValue,'vector') % vector
        vector_node = root_node.createElement('vector');
        %check the vector is row or column and set the attribute
        vector_size = size(curr_value);
        if vector_size(1) > vector_size(2) % it's row
            vector_type = 'column';
        else
            vector_type = 'row';
        end
        vector_type_attr = root_node.createAttribute('type');
        vector_type_attr.setNodeValue(vector_type);
        vector_node.setAttributeNode(vector_type_attr);
        %set the value
        for ii=1:length(curr_value)
            tmp_value = curr_value(ii);
            if isnumeric(tmp_value)
                cn_node = root_node.createElement('cn');
            else
                cn_node = root_node.createElement('ci');
            end
            cn_text = root_node.createTextNode(num2str(tmp_value));
            cn_node.appendChild(cn_text);
            vector_node.appendChild(cn_node);
        end
        param_node.appendChild(vector_node);
    elseif strcmp(class_attr.getNodeValue,'matrix') % matrix
        array_node = root_node.createElement('matrix');
        for rr = 1 : size(curr_value,1) % store the matrix row by row
            arrayrow_node = root_node.createElement('arrayrow');
            for cc = 1 : size(curr_value,2)
                tmp_value = curr_value(rr,cc);
                if isnumeric(tmp_value)
                    cn_node = root_node.createElement('cn');
                else
                    cn_node = root_node.createElement('ci');
                end
                cn_text = root_node.createTextNode(num2str(tmp_value));
                cn_node.appendChild(cn_text);
                arrayrow_node.appendChild(cn_node);
            end
            array_node.appendChild(arrayrow_node);
        end
        param_node.appendChild(array_node);
    end
    lop_node.appendChild(param_node);
end
end
