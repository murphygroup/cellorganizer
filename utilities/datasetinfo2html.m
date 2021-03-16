function answer = datasetinfo2html( model, fileID )
% DATASETINFO2HTML Generates an HTML table of the model dataset info, i.e.
% all image filenames used in preprocessing, the corresponding param files
% generated on their behalf and whether the image file was accepted for the
% purpose of generating the model itself.
%
% model     Generative Model Structure (REQUIRED)

% Nicole Matamala
%
% Copyright (C) 2019 Murphy Lab
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

answer = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK MODEL FIELDS %%%%%%%%%%%%%%%%%%%%%%%%%%%

% model.name
if ~isfield(model, 'name')
    warning("Model does not have name field. Assigning generic name 'model'.");
    name = 'model';
else
    name = string(model.name);
end

% model.dataset
if ~isfield(model, 'dataset')
    error("Model does not have dataset field.");
    return;
    % dataset.nuclear_membrane_images or dataset.cell_membrane_images
elseif ~isfield(model.dataset, 'cell_membrane_images') && ...
        ~isfield(model.dataset, 'nuclear_membrane_images')
    error("Model does not have dataset.cell_membrane image or dataset.nuclear_membrane_images field.");
    return;
    % dataset.parameterization
elseif ~isfield(model.dataset, 'parameterization') || ...
        isempty(model.dataset.parameterization)
    error("Model does not have dataset.parameterization field.");
    return;
end

%%%%%%%%%%%%%%%%%%%%%% GENERATE DATASET INFO MATRIX %%%%%%%%%%%%%%%%%%%%%%%

% Determine which columns will be built in the HTML file
if (~isfield( model.dataset, 'nuclear_membrane_images' ) || ...
        isempty( model.dataset.nuclear_membrane_images ))
    if isempty( model.dataset.cell_membrane_images )
        error("No images were in either dataset field.")
        return;
    else
        % COL 1: Cell membrane images
        cellImgs = model.dataset.cell_membrane_images;
        L = length(cellImgs);
        
        % COL 2: Sort the parameterization filenames
        params = string(sort_nat(model.dataset.parameterization));
        
        % COL 3: Status of image usage in model generation
        status = strings([1, L]);
        for i = 1:L
            [~, ~, ext] = fileparts(params(i));
            % Image file was used in the generation of the model
            if ext == '.mat'
                status(i) = 'Accepted';
                % Image file was not used in the generation of the model
            else
                status(i) = 'Rejected';
            end
        end
        
        % Tranpose matrix to have an N by 3 matrix
        MAT = transpose([cellImgs; params; status]);
        
        colName = 'Cell Membrane Files';
    end
else
    if (~isfield( model.dataset, 'cell_membrane_images' ) || ...
            isempty( model.dataset.cell_membrane_images ))
        % COL 1: Nuclear membrane images
        nucImgs = model.dataset.nuclear_membrane_images;
        L = length(nucImgs);
        
        % COL 2: Sort the parameterization filenames
        params = string(sort_nat(model.dataset.parameterization));
        
        % COL 3: Status of image usage in model generation
        status = strings([1, L]);
        for i = 1:L
            [~, ~, ext] = fileparts(params(i));
            % Image file was used in the generation of the model
            if ext == '.mat'
                status(i) = 'Accepted';
                % Image file was not used in the generation of the model
            else
                status(i) = 'Rejected';
            end
        end
        
        % Tranpose matrix to have an N by 3 matrix
        MAT = transpose([nucImgs; params; status]);
        
        colName = 'Nuclear Membrane Files';
        
    else
        % COL 1: Cell membrane images
        cellImgs = model.dataset.cell_membrane_images;
        
        % COL 2: Nuclear membrane images
        nucImgs = model.dataset.nuclear_membrane_images;
        L = length(cellImgs);
        
        % COL 3: Sort the parameterization filenames
        params = string(sort_nat(model.dataset.parameterization));
        
        % COL 4: Status of image usage in model generation
        status = strings([1, L]);
        for i = 1:L
            [~, ~, ext] = fileparts(params(i));
            % Image file was used in the generation of the model
            if ext == '.mat'
                status(i) = 'Accepted';
                % Image file was not used in the generation of the model
            else
                status(i) = 'Rejected';
            end
        end
        
        % Tranpose matrix to have an N by 4 matrix
        try
            MAT = transpose([cellImgs; nucImgs; params; status]);
        catch
            warning('Dataset is made of function handles. Ignoring dataset.');
            MAT = [];
        end
        
        colName = 'Both';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATE HTML TABLE %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initilialize HTML file
header2html(fileID, 'Dataset Information');

fprintf(fileID,['<body style= font-family:arial;>\n']);
if ~isempty(MAT)
    fprintf(fileID, ['<table border="1"; '...
        'bordercolor="#dddddd"; '...
        'cellpadding="8"; '...
        'cellspacing="0">\n']);
    if colName ~= 'Both'
        % Write HTML Table column names
        fprintf(fileID, ['<tr>\n', ...
            '<th>%s</th>\n',...
            '<th>Parameter</th>\n',...
            '<th>Preprocessing Status</th>\n',...
            '</tr>\n'], colName);
        % Write dataset info matrix to HTML table
        for i = 1:L
            fprintf(fileID, ['<tr>\n', ...
                '<td>%s</td>\n', ...
                '<td>%s</td>\n'], MAT(i, 1), MAT(i, 2));
            if MAT(i, 3) == 'Accepted'
                fprintf(fileID, '<td><font color="green">%s</font></td>\n%% </tr>\n', MAT(i, 3));
            else
                fprintf(fileID, '<td><font color="red">%s</font></td>\n%% </tr>\n', MAT(i, 3));
            end
        end
    else
        % Write HTML Table column names
        fprintf(fileID, ['<tr>\n', ...
            '<th>Cell Membrane File</th>\n',...
            '<th>Nuclear Membrane File</th>\n',...
            '<th>Parameter</th>\n',...
            '<th>Preprocessing Status</th>\n',...
            '</tr>\n']);
        % Write dataset info matrix to HTML table
        for i = 1:L
            fprintf(fileID, ['<tr>\n', ...
                '<td>%s</td>\n', ...
                '<td>%s</td>\n', ...
                '<td>%s</td>\n'], MAT(i, 1), MAT(i, 2), MAT(i, 3));
            if MAT(i, 4) == 'Accepted'
                fprintf(fileID, '<td><font color="green">%s</font></td>\n%% </tr>\n', MAT(i, 4));
            else
                fprintf(fileID, '<td><font color="red">%s</font></td>\n%% </tr>\n', MAT(i, 4));
            end
        end
    end
    
    % Close HTML table
    fprintf(fileID, ['</table>\n']);
else
    fprintf(fileID, ['Dataset is composed of function handles. Ignoring request.']);
end

answer = true;
end%datasetinfo2html