function copydemo(demo_name)
% The function is used for copy demo to a new file so that users can change
% the demos to adpat to their own work. When run the demo, it will first
% check if the demo_name is a valid demo name in CellOrganizer. Then the
% user will be asked to put in a filename for the new script. The filename
% should be a valid matlab script name. The extension ".m" is optional. And
% the will be created in the current directory. 
% 
% Author: Xiongtao Ruan xruan@andrew.cmu.edu
% Date: Jan. 7, 2016


cur_path = which(mfilename);

cellorganizer_root = fileparts(cur_path);

if nargin < 1
    error('Please provide a valid demo filename!');
end

if ~exist(demo_name, 'file') 
    error('Sorry, the demo %s  does not exist! Please check again!', demo_name);
end

% check if the demo_name is a valid name for a demo

demo_fullpath = which(demo_name);
[demo_path,demo_name_1, ext] = fileparts(demo_fullpath); 

regex_str = [cellorganizer_root, filesep, 'demos/[23]D/', demo_name_1, filesep, demo_name_1, '.m'];

if ~regexp(demo_fullpath, regex_str)
    error('It is not a valid demo! Please check again!')
end

% input a new filename 
prompt = 'What name do you want for your script?\n';
input_filename = input(prompt,'s');

new_filename = '';
while isempty(new_filename)
    if ~isempty(regexp(input_filename, '[/\*:?"<>|]', 'once'))
        prompt = sprintf('%s is not a valid filename, please try again!\n', input_filename);
        input_filename = input(prompt,'s');
    else
        if strcmp(input_filename(end - 1 : end), '.m')
            new_filename = input_filename;
        else
            new_filename = [input_filename, '.m'];
        end
    end
end

copyfile(demo_fullpath, new_filename);

% edit(new_filename);

end


