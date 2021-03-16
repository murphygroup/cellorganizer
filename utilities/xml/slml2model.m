function list_cell = slml2model( filename )
%SLML2MODEL Reads an SLML Level 1 instance and returns a valid model that can be used in CellOrganizer.
%
%Inputs
%filename: a valid SLML instance
%
%Outputs
%model: a Matlab structure holding the generative model described in the SLML instance
%
%See also MODEL2SLML

% Author: Yue Yu
% 7/11/2012 Modified by Y.Yu
%
% Copyright (C) 2012-2014 Murphy Lab
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
%
% 7/16/2012 Y.Yu Fixed a bug when reading 2D model, now it can read models
% which contain cell array


%step 1: check input and output arguments
if nargin ~= 1
    error('Wrong number of input arguments.')
end

model = struct;
if ~exist( filename, 'file' )
    warning( ['Input argument ' filename ' does not exist.'] );
    return
end

%step 2: read slml instance to DOM
DOMnode = xmlread(filename);

%step 3: parse to matlab structure
docRoot = DOMnode.getDocumentElement(); %root node of DOM

java_listOfCells = docRoot.getElementsByTagName('ListOfCells').item(0);
java_list_cell = java_listOfCells.getElementsByTagName('cell');
numberOfCells = java_list_cell.getLength;
list_cell = {};

model = struct;
for cell_i=0:numberOfCells-1
    java_cell = java_list_cell.item(cell_i);
    
    java_listOfModels = java_cell.getElementsByTagName('ListOfModels').item(0);
    java_list_model = java_listOfModels.getElementsByTagName('model');
    numberOfModels = java_list_model.getLength;
    list_model = {};
    for model_i=0:numberOfModels-1
        java_model = java_list_model.item(model_i);
        
        java_listOfCompartments = java_model.getElementsByTagName('ListOfCompartment').item(0);
        java_list_compartment = java_listOfCompartments.getElementsByTagName('compartment');
        numberOfCompartments = java_list_compartment.getLength;
        list_compartment = [];
        %set model attributes
        model_attr = java_model.getAttributes;
        model_attr_num = model_attr.getLength;
        for model_attr_i = 0 : model_attr_num-1
            model.(model_attr.item(model_attr_i).getNodeName.toCharArray') = ...
            model_attr.item(model_attr_i).getNodeValue.toCharArray';
        end
        
        for compartment_i=0:numberOfCompartments-1
            java_compartment = java_list_compartment.item(compartment_i);
            %get all compartment attributes
            comp_attr = java_compartment.getAttributes;
            comp_attr_num = comp_attr.getLength;
            compartment_attributes = {};% store as {name,value,name,value...}
            for comp_attr_i = 0 : comp_attr_num-1
                compartment_attributes{end+1} = comp_attr.item(comp_attr_i).getNodeName.toCharArray';
                compartment_attributes{end+1} = comp_attr.item(comp_attr_i).getNodeValue.toCharArray';
            end
            comp_attr_name_idx = find(ismember(compartment_attributes,'name')>0);
            comp_attr_name = compartment_attributes{comp_attr_name_idx+1};
            %remove the name attr from the cell array
            compartment_attributes(comp_attr_name_idx) = [];
            compartment_attributes(comp_attr_name_idx) = [];
            
            % check if it's protein model, nuclear model or cell model
            if ~isempty(strfind(comp_attr_name,'nuclear')) % it's nuclear shape model
                java_shape = java_compartment.getElementsByTagName('shape').item(0);          
                java_listOfParameters = java_shape.getElementsByTagName('ListOfParameters').item(0);
                java_list_parameter = java_listOfParameters.getElementsByTagName('parameter');
                numberOfParameters = java_list_parameter.getLength;
                list_parameter = {};
                for parameter_i=0:numberOfParameters-1
                   java_parameter = java_list_parameter.item(parameter_i);
                   [parameter,param_name] = makeParameter(java_parameter);
                   list_parameter{end+1} = param_name;
                   list_parameter{end+1} = parameter;
                end
                curr_lop = mergeparameter(list_parameter);% merge the parameters, make the lop
                %set all lop attributes as parameter
                attr_param_L = java_listOfParameters.getAttributes.getLength; % the number of attr parameters
                for attr_param_i = 0 : attr_param_L-1
                    curr_lop.(java_listOfParameters.getAttributes.item(attr_param_i).getNodeName.toCharArray') = ...
                    java_listOfParameters.getAttributes.item(attr_param_i).getNodeValue.toCharArray';
                end
                
                model.(comp_attr_name) = curr_lop;
                for i = 1 : length(compartment_attributes)-1
                    model.(comp_attr_name).(compartment_attributes{i}) = compartment_attributes{i+1};
                end
            elseif ~isempty(strfind(comp_attr_name,'cell')) % it's cell shape model
                java_shape = java_compartment.getElementsByTagName('shape').item(0);
                
                java_listOfParameters = java_shape.getElementsByTagName('ListOfParameters').item(0);
                java_list_parameter = java_listOfParameters.getElementsByTagName('parameter');
                numberOfParameters = java_list_parameter.getLength;
                list_parameter = {};
                for parameter_i=0:numberOfParameters-1
                   java_parameter = java_list_parameter.item(parameter_i);
                   [parameter,param_name] = makeParameter(java_parameter);
                   list_parameter{end+1} = param_name;
                   list_parameter{end+1} = parameter;
                end
                curr_lop = mergeparameter(list_parameter);% merge the parameters, make the lop
                %set all lop attributes as parameter
                attr_param_L = java_listOfParameters.getAttributes.getLength; % the number of attr parameters
                for attr_param_i = 0 : attr_param_L-1
                    curr_lop.(java_listOfParameters.getAttributes.item(attr_param_i).getNodeName.toCharArray') = ...
                    java_listOfParameters.getAttributes.item(attr_param_i).getNodeValue.toCharArray';
                end
                
                model.(comp_attr_name) = curr_lop;
                for i = 1 : length(compartment_attributes)-1
                    model.(comp_attr_name).(compartment_attributes{i}) = compartment_attributes{i+1};
                end
            else %it's protein model
                if java_compartment.getElementsByTagName('shape').getLength > 0
                    java_shape = java_compartment.getElementsByTagName('shape').item(0);
                
                    java_listOfParameters = java_shape.getElementsByTagName('ListOfParameters').item(0);
                    java_list_parameter = java_listOfParameters.getElementsByTagName('parameter');
                    numberOfParameters = java_list_parameter.getLength;
                    list_parameter = {};
                    for parameter_i=0:numberOfParameters-1
                       java_parameter = java_list_parameter.item(parameter_i);
                       [parameter,param_name] = makeParameter(java_parameter);
                       list_parameter{end+1} = param_name;
                       list_parameter{end+1} = parameter;
                    end
                    curr_lop = mergeparameter(list_parameter);% merge the parameters, make the lop
                    %set all lop attributes as parameter
                    attr_param_L = java_listOfParameters.getAttributes.getLength; % the number of attr parameters
                    for attr_param_i = 0 : attr_param_L-1
                        curr_lop.(java_listOfParameters.getAttributes.item(attr_param_i).getNodeName.toCharArray') = ...
                        java_listOfParameters.getAttributes.item(attr_param_i).getNodeValue.toCharArray';
                    end
                    if strcmp(model.dimensionality,'2D')
                        model.(comp_attr_name).objectModel = curr_lop;
                    else
                        model.(comp_attr_name).size = curr_lop;
                    end
                end
                if java_compartment.getElementsByTagName('texture').getLength > 0
                    java_texture = java_compartment.getElementsByTagName('texture').item(0);
                
                    java_listOfParameters = java_texture.getElementsByTagName('ListOfParameters').item(0);
                    java_list_parameter = java_listOfParameters.getElementsByTagName('parameter');
                    numberOfParameters = java_list_parameter.getLength;
                    list_parameter = {};
                    for parameter_i=0:numberOfParameters-1
                       java_parameter = java_list_parameter.item(parameter_i);
                       [parameter,param_name] = makeParameter(java_parameter);
                       list_parameter{end+1} = param_name;
                       list_parameter{end+1} = parameter;
                    end
                    curr_lop = mergeparameter(list_parameter);% merge the parameters, make the lop
                    %set all lop attributes as parameter
                    attr_param_L = java_listOfParameters.getAttributes.getLength; % the number of attr parameters
                    for attr_param_i = 0 : attr_param_L-1
                        curr_lop.(java_listOfParameters.getAttributes.item(attr_param_i).getNodeName.toCharArray') = ...
                        java_listOfParameters.getAttributes.item(attr_param_i).getNodeValue.toCharArray';
                    end
                
                    model.(comp_attr_name).texture = curr_lop;
                end
                if java_compartment.getElementsByTagName('position').getLength > 0
                    java_position = java_compartment.getElementsByTagName('position').item(0);
                
                    java_listOfParameters = java_position.getElementsByTagName('ListOfParameters').item(0);
                    java_list_parameter = java_listOfParameters.getElementsByTagName('parameter');
                    numberOfParameters = java_list_parameter.getLength;
                    list_parameter = {};
                    for parameter_i=0:numberOfParameters-1
                       java_parameter = java_list_parameter.item(parameter_i);
                       [parameter,param_name] = makeParameter(java_parameter);
                       list_parameter{end+1} = param_name;
                       list_parameter{end+1} = parameter;
                    end
                    curr_lop = mergeparameter(list_parameter);% merge the parameters, make the lop
                    %set all lop attributes as parameter
                    attr_param_L = java_listOfParameters.getAttributes.getLength; % the number of attr parameters
                    for attr_param_i = 0 : attr_param_L-1
                        curr_lop.(java_listOfParameters.getAttributes.item(attr_param_i).getNodeName.toCharArray') = ...
                        java_listOfParameters.getAttributes.item(attr_param_i).getNodeValue.toCharArray';
                    end
                    if strcmp(model.dimensionality,'2D')
                        model.(comp_attr_name).positionModel = curr_lop;
                    else
                        model.(comp_attr_name).position = curr_lop;
                    end
                end
                if java_compartment.getElementsByTagName('frequency').getLength > 0
                    java_frequency = java_compartment.getElementsByTagName('frequency').item(0);
                
                    java_listOfParameters = java_frequency.getElementsByTagName('ListOfParameters').item(0);
                    java_list_parameter = java_listOfParameters.getElementsByTagName('parameter');
                    numberOfParameters = java_list_parameter.getLength;
                    list_parameter = {};
                    for parameter_i=0:numberOfParameters-1
                       java_parameter = java_list_parameter.item(parameter_i);
                       [parameter,param_name] = makeParameter(java_parameter);
                       list_parameter{end+1} = param_name;
                       list_parameter{end+1} = parameter;
                    end
                    curr_lop = mergeparameter(list_parameter);% merge the parameters, make the lop
                    %set all lop attributes as parameter
                    attr_param_L = java_listOfParameters.getAttributes.getLength; % the number of attr parameters
                    for attr_param_i = 0 : attr_param_L-1
                        curr_lop.(java_listOfParameters.getAttributes.item(attr_param_i).getNodeName.toCharArray') = ...
                        java_listOfParameters.getAttributes.item(attr_param_i).getNodeValue.toCharArray';
                    end
                    
                    model.(comp_attr_name).frequency = curr_lop;
                end
                for i = 0 : length(compartment_attributes)/2-1
                    model.(comp_attr_name).(compartment_attributes{i*2+1}) = compartment_attributes{i*2+2};
                end
            end     
        end
        list_model{end+1} = model;
        
    end
    list_cell{end+1} = list_model;
end

%step 5: enjoy!
end


%%helper functions
function lop = mergeparameter(rawlist)
% merge all parameters with the same field names
% @param rawlist, cell array, {topname,value,topname,value...}, each value is a struct
% @return a struct, fields with unique fields names from all topnames

N = length(rawlist)/2;% the total number of parameters
if N ~= floor(N)
    error ('wrong parameter list');
end
namelist = {};
for ii = 0 : N-1
    namelist{end+1} = rawlist{ii*2+1};
end
uni_names = unique(namelist);% the unique name list
uni_num = length(uni_names);% the number of unique names
for ii = 1 : uni_num
    tmp_name = uni_names{ii};
    tmp_idx = find(ismember(namelist,tmp_name)>0); % the index of current field name
    if length(tmp_idx) == 1 % this field doesn't need to be merged
        lop.(tmp_name) = rawlist{tmp_idx(1)*2};
    else % this field needs to be merged
        
        %first check if their fields have save name
        subfield_names = {};
        subfield_value = {};
        for jj = 1 : length(tmp_idx)
            tmp_value = rawlist{tmp_idx(jj)*2};
            subname = fieldnames(tmp_value);
            subfield_names{end+1} = subname{1}; % store all sub field names
            subfield_value{end+1} = tmp_value; % store all same top name values
        end
        
        if length(unique(subfield_names)) == length(subfield_names)% all sub fields are differnet
            for jj = 1 : length(tmp_idx)
                lop.(tmp_name).(subfield_names{jj}) = getfield(rawlist{tmp_idx(jj)*2},subfield_names{jj});
            end    
        else % some sub field also need to be merged, use the recursion!
            sub_unique = unique(subfield_names);
            num_subuni = length(sub_unique);% the number of unique sub fields
            % merge all same sub field and return as one sub field
            for jj = 1 : num_subuni
                tmp_subname = sub_unique{jj};
                tmp_subidx = find(ismember(subfield_names,tmp_subname)>0); % the index of current field name
                if length(tmp_subidx) >= 2 %it needs to be merged
                    tmp_rawlist = {}; % make the rawlist for recursion
                    for k = 1 : length(tmp_subidx)
                        tmp_rawlist{end+1} = tmp_subname;
                        tmp_rawlist{end+1} = getfield(subfield_value{tmp_subidx(k)},tmp_subname);
                    end
                    tmp_lop = mergeparameter(tmp_rawlist);% recursion
                    lop.(tmp_name).(tmp_subname) = getfield(tmp_lop,tmp_subname);
                else % it doesn't have to be merged
                    lop.(tmp_name).(tmp_subname) = getfield(subfield_value{tmp_subidx(1)},tmp_subname);
                end
            end
        end
        %convert the struct back to cell array        
        iscellstruct = 1;% mark if it's cellstruct
        cellL = 0;
        struct2cell_L = length(subfield_names);
        for cellstruct_i = 1 : struct2cell_L
            if length(subfield_names{cellstruct_i}) <= 10 
                iscellstruct = 0;
            else
                if ~strcmp(subfield_names{cellstruct_i}(1:10),'cellstruct')
                    iscellstruct = 0;
                end
            end 
            if str2double(subfield_names{cellstruct_i}(11:end)) > cellL
                cellL = str2double(subfield_names{cellstruct_i}(11:end));
            end
        end
        tmp_cell = {};
        if iscellstruct % the tmp_name field need to be converted
            for cellstruct_i = 1 : cellL
                tmp_cell{cellstruct_i} = getfield(lop.(tmp_name),strcat('cellstruct',num2str(cellstruct_i)));
            end
            lop = rmfield(lop,tmp_name);
            lop.(tmp_name) = tmp_cell;
        end       
    end
end
end

function [parameter,topname] = makeParameter(java_parameter)
%  Converts a parameter DOM object into a struct, the structure of the struct
%  is determined by the parameter's name
%  @param the parameter DOM object(java_parameter)
%  @return a matlab struct(parameter), and the variable name of the struct
     numofattr = java_parameter.getAttributes.getLength;
     attr_names = {}; % used for storing the attr names
     attr_values = {}; % used for storing the attr values
     % retrieve all attributes
     for i = 0 : numofattr-1
         attr_names{end+1} = java_parameter.getAttributes.item(i).getNodeName.toCharArray';
         attr_values{end+1} = java_parameter.getAttributes.item(i).getNodeValue.toCharArray';
     end
     % find the class attribute
     [~,class_idx] = ismember('class',attr_names);
     param_class = attr_values{class_idx}; % get the parameter's class
     % find the isNumeric attribute
     [~,numeric_idx] = ismember('isNumeric',attr_names);
     param_numeric = attr_values{numeric_idx}; % get the parameter's isNumeric
     % find the type attribute
     [~,type_idx] = ismember('type',attr_names);
     param_type = attr_values{type_idx}; % get the parameter's type
     % find the name attribute
     [~,name_idx] = ismember('name',attr_names);
     param_name = attr_values{name_idx}; % get the parameter's name
     
     %read all parameters
     if strcmp(param_class,'scalar') % the parameter is scalar
         if strcmp(param_numeric,'true') % the parameter is not string
             param_value = str2num(java_parameter.getElementsByTagName('scalar').item(0).getTextContent.toCharArray');
         else % the parameter is a string
             param_value = java_parameter.getElementsByTagName('scalar').item(0).getTextContent.toCharArray';
         end
     elseif strcmp(param_class,'vector') % the parameter is a vector
         try
             Lofvec = ...
             java_parameter.getElementsByTagName('vector').item(0).getElementsByTagName('cn').getLength; 
             % get the length of the vector
         catch
             Lofvec = ...
             java_parameter.getElementsByTagName('vector').item(0).getElementsByTagName('ci').getLength; 
             % get the length of the vector
         end
         param_value = [];% initialize the vector
         if strcmp(param_numeric,'true') % the parameter is not string
             for i = 0 : Lofvec-1
                 param_value(end+1) = ...
                 str2num(java_parameter.getElementsByTagName('vector').item(0).getElementsByTagName('cn').item(i).getTextContent);
             end
         else % the parameter is string
             param_value = {};
             for i = 0 : Lofvec-1
                 param_value{end+1} = ...
                 java_parameter.getElementsByTagName('vector').item(0).getElementsByTagName('ci').item(i).getTextContent;
             end
         end
         try
             vector_type = java_parameter.getElementsByTagName('vector').item(0).getAttributes.item(0).getNodeValue.toCharArray';
             if strcmp(vector_type,'column')
                 param_value = param_value';
             end
         catch
         end
           
             
     elseif strcmp(param_class,'matrix') % the parameter is a matrix
         %get the number of rows and columns
         numofrows = java_parameter.getElementsByTagName('matrix').item(0).getElementsByTagName('arrayrow').getLength;
         try
             numofcol = java_parameter.getElementsByTagName('matrix').item(0).getElementsByTagName('arrayrow').item(0).getElementsByTagName('cn').getLength;
         catch
             numofcol = java_parameter.getElementsByTagName('matrix').item(0).getElementsByTagName('arrayrow').item(0).getElementsByTagName('ci').getLength;       
         end
         param_value = zeros(numofrows,numofcol);
         if strcmp(param_numeric,'true') % the parameter is not string
             for i = 0 : numofrows-1
                 for j = 0 : numofcol-1
                     tmp_row = java_parameter.getElementsByTagName('matrix').item(0).getElementsByTagName('arrayrow').item(i);
                     param_value(i+1,j+1) = str2num(tmp_row.getElementsByTagName('cn').item(j).getTextContent);
                 end
             end
         else % the parameter is string
             param_value = cell(numofrows,numofcol);
             for i = 0 : numofrows-1
                 for j = 0 : numofcol-1
                     tmp_row = java_parameter.getElementsByTagName('matrix').item(0).getElementsByTagName('arrayrow').item(i);
                     param_value{i+1,j+1} = tmp_row.getElementsByTagName('ci').item(j).getTextContent;
                 end
             end
         end
     else
         error ('Unknown class')
     end
     
     %split the name
     split_name = regexp(param_name,'\.','split');
     nameL = length(split_name); % get the length of name
     if nameL == 3
         topname = split_name{2};
         parameter.(split_name{nameL}) = param_value;
     elseif nameL == 2
         topname = split_name{2};
         parameter = param_value;
     elseif nameL >= 4
         topname = split_name{2};
         parameter.(split_name{nameL}) = param_value;
         for i = nameL-1 : -1 : 3
             rmname = fieldnames(parameter);
             parameter = setfield(parameter,split_name{i},parameter);
             parameter = rmfield(parameter,rmname);
         end
     else
         error('wrong name is presented')
     end
 end



