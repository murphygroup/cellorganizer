function answer = segmentation2html( model, fileID, save_folder )
% SEGMENTATION2HTML Generates an HTML table of the segmented cell and
% nucleus images.
%
% model     Generative Model Structure (REQUIRED)

% N. Matamala
%
% Copyright (C) 2018 Murphy Lab
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
elseif ~isfield(model.dataset, 'segmentation')
    error("Model does not have dataset.segmentation field.");
    return;
elseif isempty( model.dataset.segmentation )
    error("Model does not have segmentation images within segmentation field.");
    return;
end

%%%%%%%%%%%%%%%%%%%%%% GENERATE GALLERY IMAGE MATRIX %%%%%%%%%%%%%%%%%%%%%%

% Create directory for images
if ~exist( [ pwd filesep save_folder '/images'], 'dir' )
    mkdir( [save_folder '/images'] );
end

if isfield( model.dataset.segmentation{1, 1}, 'cell' ) && ...
   isfield( model.dataset.segmentation{1, 1}, 'nucleus' )

    cellfName = model.dataset.cell_membrane_images;
    nucfName = model.dataset.nuclear_membrane_images;
    
    L = length(nucfName);
    
    % Generate cell and nucleus images
    for i = 1:L
        if isstruct(model.dataset.segmentation{i})
            cellImg = im2projection(model.dataset.segmentation{1, i}.cell);
            nucImg = im2projection(model.dataset.segmentation{1, i}.nucleus);
            imwrite(cellImg, strcat(save_folder, '/images/cell', string(i), '.jpeg'));
            imwrite(nucImg, strcat(save_folder, '/images/nuc', string(i), '.jpeg'));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%% GENERATE HTML TABLE %%%%%%%%%%%%%%%%%%%%%%%%%

    % Initilialize HTML file
    htmlT = fileID;
    fprintf(htmlT, '%%%% Segmentation\n');
    fprintf(htmlT,['%% <html>\n']);
    fprintf(htmlT,['%% <body style= font-family:arial;>\n']);

    fprintf(htmlT, ['%% <table border="1"; '...
    'bordercolor="#dddddd"; '...
    'cellpadding="8"; '...
    'cellspacing="0">\n']);

    % Write filename and segmentation images to HTML table
    for i = 1:L
        if isstruct(model.dataset.segmentation{i})
            fprintf(htmlT, ['%% <tr>\n', ...
            '%% <th colspan="2">%s</th>\n', ...
            '%% </tr>\n'], strcat("Image ", string(i)));
            fprintf(htmlT, ['%% <tr>\n', ...
            '%% <td><center>Cell Membrane</center></td>\n', ...
            '%% <td><center>Nuclear Membrane</center></td>\n', ...
            '%% </tr>\n']);
            fprintf(htmlT, ['%% <tr>\n',...
            '%% <td>%s</td>\n', ...
            '%% <td>%s</td>\n', ...
            '%% </tr>\n'], string(cellfName(i)), string(nucfName(i)));
            cellName = strcat('images/cell', string(i), '.jpeg');
            nucName = strcat('images/nuc', string(i), '.jpeg');
            fprintf(htmlT, ['%% <tr>\n', ...
            '%% <td><center><a href="%s"><img src="%s" height="256" width="256"></a></td>\n', ...
            '%% <td><center><a href="%s"><img src="%s" height="256" width="256"></a></td>\n', ...
            '%% </tr>\n'], cellName, cellName, nucName, nucName);
        end
    end
    % Close HTML table
    fprintf(htmlT, ['%% </table>\n%% </body>\n']);

    % Close HTML file
    fprintf(htmlT, ['%% </html>\n']);
    
    answer = true;

elseif isfield( model.dataset.segmentation{1, 1}, 'cell' ) && ...
       ~isfield( model.dataset.segmentation{1, 1}, 'nucleus' )
    filenames = model.dataset.cell_membrane_images;
    L = length(filenames);
    
    % Generate cell and nucleus images
    for i = 1:L
        if isstruct(model.dataset.segmentation{i})
            cellImg = im2projection(model.dataset.segmentation{1, i}.cell);
            imwrite(cellImg, strcat(save_folder,'/images/cell', string(i), '.jpeg'));
        end
    end
    
    % Initilialize HTML file
    htmlT = fileID;
    fprintf(htmlT, '%%%% Segmentation\n');
    fprintf(htmlT,['%% <html>\n']);
    fprintf(htmlT,['%% <body style= font-family:arial;>\n']);

    fprintf(htmlT, ['%% <table border="1"; '...
    'bordercolor="#dddddd"; '...
    'cellpadding="8"; '...
    'cellspacing="0">\n']);

    % Write filename and segmentation images to HTML table
    for i = 1:L
        if isstruct(model.dataset.segmentation{i})
            fprintf(htmlT, ['%% <tr>\n', ...
            '%% <th>%s</th>\n', ...
            '%% </tr>\n'], strcat("Image ", string(i)));
            fprintf(htmlT, ['%% <tr>\n', ...
            '%% <td><center>Cell Membrane</center></td>\n', ...
            '%% </tr>\n']);
            fprintf(htmlT, ['%% <tr>\n',...
            '%% <td>%s</td>\n', ...
            '%% </tr>\n'], string(filenames(i)));
            cellName = strcat('images/cell', string(i), '.jpeg');
            fprintf(htmlT, ['%% <tr>\n', ...
            '%% <td><center><a href="%s"><img src="%s" height="256" width="256"></a></td>\n', ...
            '%% </tr>\n'], cellName, cellName);
        end
    end
    % Close HTML table
    fprintf(htmlT, ['%% </table>\n%% </body>\n']);

    % Close HTML file
    fprintf(htmlT, ['%% </html>\n']);
    
    answer = true;
    
elseif ~isfield( model.dataset.segmentation{1, 1}, 'cell' ) && ...
        isfield( model.dataset.segmentation{1, 1}, 'nucleus' )
    filenames = model.dataset.nuclear_membrane_images;
    L = length(filenames);
    
    % Generate cell and nucleus images
    for i = 1:L
        if isstruct(model.dataset.segmentation{i})
            nucImg = im2projection(model.dataset.segmentation{1, i}.nucleus);
            imwrite(nucImg, strcat(save_folder, '/images/nuc', string(i), '.jpeg'));
        end
    end
    
    % Initilialize HTML file
    htmlT = fileID;
    fprintf(htmlT, '%%%% Segmentation\n');
    fprintf(htmlT,['%% <html>\n']);
    fprintf(htmlT,['%% <body style= font-family:arial;>\n']);

    fprintf(htmlT, ['%% <table border="1"; '...
    'bordercolor="#dddddd"; '...
    'cellpadding="8"; '...
    'cellspacing="0">\n']);

    % Write filename and segmentation images to HTML table
    for i = 1:L
        if isstruct(model.dataset.segmentation{i})
            fprintf(htmlT, ['%% <tr>\n', ...
            '%% <th>%s</th>\n', ...
            '%% </tr>\n'], strcat("Image ", string(i)));
            fprintf(htmlT, ['%% <tr>\n', ...
            '%% <td><center>Nuclear Membrane</center></td>\n', ...
            '%% </tr>\n']);
            fprintf(htmlT, ['%% <tr>\n',...
            '%% <td>%s</td>\n', ...
            '%% </tr>\n'], string(filenames(i)));
            nucName = strcat('images/nuc', string(i), '.jpeg');
            fprintf(htmlT, ['%% <tr>\n', ...
            '%% <td><center><a href="%s"><img src="%s" height="256" width="256"></a></td>\n', ...
            '%% </tr>\n'], nucName, nucName);
        end
    end
    % Close HTML table
    fprintf(htmlT, ['%% </table>\n%% </body>\n']);

    % Close HTML file
    fprintf(htmlT, ['%% </html>\n']);
    
    answer = true;
end

end