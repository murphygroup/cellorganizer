function [result] = get_synapse_file_data(filename, include_cell_index)
% Returns a cell array where index is timeframe, and each element contains an array where a row is a cell (not necessarily the same one at the same index across timeframes) and columns are X, Y, and T (timeframe relative to synapse formation (T=0)).
% Just read the whole synapse file and return its parsed contents.
% 2012-09-02 tebuck: Copied from BK's read_synapse_file_new127.m.
  
  % error('Implementation yet unfinished below this line!')
  
if ~exist('include_cell_index', 'var')
  include_cell_index = false;
end

% fname = filename{1};
fname = filename;
%fname, pause
fid=fopen(fname);

slash_ = find(fname == '/');

Cell_index = 0;
Flag_new_cell = false;

cell_fname = '';
final_exists = false;

synapses_per_timeframe = {};

while ~feof(fid)
    tline = fgetl(fid);
    if isempty(tline)
        continue;
    end
    if tline(1) == '#'  % ignore this line & recognize a new cell starts from the next line
        Cell_index = Cell_index + 1;
        Flag_new_cell = true;
        Cell_volume_first = 0;
        Cell_fluo_first = 0;
        
    else
        [T, X, Y, timeframe, Manual, opt_seg_1, opt_seg_2, opt_cutoff, opt_fg_seed, opt_fg_seed2]=strread(tline,'%d %d %d %d %c %f %d %f %d %d');
        if isempty(X) || isempty(Y) || isempty(timeframe)
            continue;
        else
            % Check if this work has been done:
            
            while length(synapses_per_timeframe) < timeframe
              if include_cell_index
                synapses_per_timeframe{end + 1} = zeros(0, 4);
              else
                synapses_per_timeframe{end + 1} = zeros(0, 3);
              end
            end
            
            % The index into synapses_per_timeframe gives timeframe:
            if include_cell_index
              synapses_per_timeframe{timeframe} = [synapses_per_timeframe{timeframe}; ...
                X, Y, T, Cell_index];
            else
              synapses_per_timeframe{timeframe} = [synapses_per_timeframe{timeframe}; ...
                X, Y, T];
            end
            % synapses_per_timeframe{timeframe} = [synapses_per_timeframe{timeframe}; ...
              % X, Y, timeframe];

          % st=['T: ' num2str(T) ' X: ' num2str(X) ' Y: ' num2str(Y) ' t: ' num2str(timeframe)];
          % disp(st);
        end
    end
end


fclose('all');


result = synapses_per_timeframe;

