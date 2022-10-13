function is_deployed( varargin )

disp('Checking number of input arguments')
if length(varargin) == 1
    text_file = varargin{1};
else
    error('Deployed function takes only 1 argument. Exiting method.');
    return
end

disp('Checking existence of input file')
[filepath, name, ext] = fileparts(text_file);

if ~exist(text_file, 'file')
    warning('Input file does not exist. Exiting method.');
    return
end

disp(['Attempting to read input file ' text_file]);
fid = fopen(text_file, 'r' );

disp('Evaluating lines from input file');
while ~feof(fid)
    line = fgets(fid);
    disp(line);
    try
        eval(line);
    catch err
        disp('Unable to parse line');
        getReport(err)
        return
    end

disp('Closing input file')
fclose(fid);

if ~exist('options', 'var')
    options = {};
end

end