function chunk_finish(path, fname)
%[can_start, final_name, final_exists] = chunk_start(path, fname, final_extension, should_block)
%chunk_finish(path, fname, unique_filename)
%  chunk_start checks if a certain unit of work has been done or is being 
%  computed. The path is the location of both the final and temporary/lock
%  files. Final files are assumed to be named [fname, '.mat'] by default,
%  but the extension can be changed with the optional argument
%  final_extension. Temporary files are named [fname, '.tmp']. Save to
%  final_name and then call chunk_finish to delete the temporary file.
%
%  Copyright 2008-2013 Taraz Buck/tebuck at cmu.
%
%   Aug 2, 2013 - G. Johnson - Added flexability for path names  


  
if exist('fname', 'var')
    temp_name = [path, '/', fname];
else
  temp_name = path;
end

if ~strcmpi(temp_name(end-3:end), '.tmp')
    temp_name = [temp_name '.tmp'];
end

if exist(temp_name, 'file')
    delete(temp_name);
end

