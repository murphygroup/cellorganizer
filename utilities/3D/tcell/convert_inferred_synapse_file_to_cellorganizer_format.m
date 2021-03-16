function [reformat_csv_filename] = convert_inferred_synapse_file_to_cellorganizer_format(inferred_csv_filename, run_folder, savename)
% convert inferred synapse file to the original format so that the pipeline
% can read it. 

if nargin < 1
    inferred_csv_filename = './inferred_synapse_location/Result_22.csv';
    run_folder = '../../../images/LAT/5C.C7 LAT 15 03 16/run 1/';
    savename = './inferred_synapse_location/test.csv';    
end

% [~, filename] = fileparts(inferred_csv_filename);
% reformat_csv_filename = sprintf('%s.csv', filename);
reformat_csv_filename = savename;

dir_info = dir(run_folder);
dir_info(1 : 2) = [];
dir_name = {dir_info.name};
is_gfp_dirname = regexpi(dir_name, '.*gfp.*');
gfp_dirname = dir_name(~cellfun(@isempty, is_gfp_dirname));
gfp_dirname = gfp_dirname{1};
gfp_full_dirname = [run_folder, '/', gfp_dirname, '/'];
gfp_full_dirname = clean_absolute_path(gfp_full_dirname);
gfp_dir_info = dir(gfp_full_dirname);
gfp_dir_info(1 : 2) = [];
gfp_filenames = {gfp_dir_info.name};

run_data_mat = csvread(inferred_csv_filename, 1, 0);

cell_ind = 0;
run_data_cell = {};
run_raw_data_cell = {};
new_cell_flag = true;

for i = 1 : size(run_data_mat, 1)
    cur_line = run_data_mat(i, :);
    if all(cur_line == -1) 
        if ~new_cell_flag 
            new_cell_flag = true;
        end
        continue;
    end    
    if new_cell_flag
        cell_ind = cell_ind + 1;
        run_data_cell{cell_ind} = {};
        run_raw_data_cell{cell_ind} = {};
        new_cell_flag = false;
    end
    if all(cur_line(1:4) == 0)
        continue;
    end
    
    cur_reformat_line = {}; 
    
    Lx = cur_line(1);
    Ly = cur_line(2);
    Rx = cur_line(3);
    Ry = cur_line(4);
    cur_frame = cur_line(5);
    cur_relative_time = cur_line(5) - cur_line(6);
    
    % [Distance, Angle, Left, Top, Width, Height] = convert_coordinates_function(Lx, Ly, Rx, Ry);
    
    cur_reformat_line{1} = [gfp_full_dirname, gfp_filenames{cur_frame}];
    cur_reformat_line{2} = 1;
    cur_reformat_line{3} = Lx;
    cur_reformat_line{4} = Ly;
    cur_reformat_line{5} = Rx;
    cur_reformat_line{6} = Ry;
    cur_reformat_line{7} = cur_relative_time;
    run_data_cell{cell_ind}{end+1,1} = cur_reformat_line;
    run_raw_data_cell{cell_ind}{end+1, 1} = cur_line;   
    
    if false
        position_info = [Distance, Angle, Left, Top, Width, Height];
        [L_x, L_y, R_x, R_y] = convert_position_to_coordinates(position_info);
        cur_line
        restored_line = [L_x, L_y, R_x, R_y]
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write to csv file

fid = fopen(reformat_csv_filename, 'wt');
header_string = 'Filename, Channel, L_X, L_Y, R_X, R_Y, Relative time';

empty_cell = repmat({''}, 1, numel(run_data_cell{1}{1}));

fwrite(fid, header_string);
fprintf(fid, '\n');

for i = 1 : numel(run_data_cell)
    cur_cell_data = run_data_cell{i};
    
    if i > 1
        for k = 1 : 2
            fwrite(fid, strjoin(empty_cell, ','));
            fprintf(fid, '\n');
        end
    end 
    
    for j = 1 : numel(cur_cell_data)
        cur_line = cur_cell_data{j};         
        for k = 1 : numel(cur_line)
            if k > 1
                fprintf(fid, ',');
            end
            if ischar(cur_line{k})
                fwrite(fid, cur_line{k});    
            else
                fprintf(fid, '%f', cur_line{k});
            end
        end
        fprintf(fid, '\n');
    end
end

fclose(fid);

end


function [Distance, Angle, Left, Top, Width, Height] = convert_coordinates_function(Lx, Ly, Rx, Ry)

    Distance = sqrt((Lx - Rx) .^ 2 + (Ly - Ry) .^ 2); 
    Width = abs(Rx - Lx);
    Height = abs(Ry - Ly);
    x = (Lx + Rx) / 2;
    y = (Ly + Ry) / 2;
    Left = x - Width / 2;
    Top = y - Height / 2;
    offset_1 = x - Lx;
    offset_2 = y - Ly;    
    Angle = -atan2(offset_2, offset_1) * 180 / pi;

end


function [L_x, L_y, R_x, R_y] = convert_position_to_coordinates(position_info)

Distance = position_info(1);
Angle = position_info(2);
Left = position_info(3);
Top = position_info(4);
Width = position_info(5);
Height = position_info(6);

x = Left + Width / 2;
y = Top + Height / 2;
seg_offset = [cos(deg2rad(-Angle)), sin(deg2rad(-Angle))] * Distance * 0.5;

L_x = x - seg_offset(1);
L_y = y - seg_offset(2);
R_x = x + seg_offset(1);
R_y = y + seg_offset(2);

end
