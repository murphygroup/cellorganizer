function [can_start, final_name, final_exists, temp_exists] = chunk_start_clean(path, fname, final_extension, should_block)
  %[can_start, final_name, final_exists] = chunk_start_clean(path, fname, final_extension, should_block)
  %chunk_finish(path, fname, unique_filename)
  %  chunk_start checks if a certain unit of work has been done or is being 
  %  computed. The path is the location of both the final and temporary/lock
  %  files. Final files are assumed to be named [fname, '.mat'] by default,
  %  but the extension can be changed with the optional argument
  %  final_extension. Temporary files are named [fname, '.tmp']. Save to
  %  final_name and then call chunk_finish to delete the temporary file.
  %
  %
  %
  % 2016-02-24 xruan: copied from chunk_start.m
  % the purpose of this file is to make a clean mechanism, we store the job
  % id in the tem file, if the temp file exist but the job isn't running,
  % then assign the job for this job. Also simplify the temp file creat
  % procedure. That is directly make the file. 
  
  
  % Process options
  if nargin < 3 || isempty(final_extension)
      final_extension = '.mat';
  end
  if nargin < 4 || isempty(should_block)
      should_block = false;
  end
  temp_name = [path, '/', fname, '.tmp'];
  final_name = [path, '/', fname, final_extension];
  
  if ~exist(temp_name, 'file') && exist(final_name, 'file')
    can_start = false;
    final_exists = true;
    temp_exists = false;
    return;
    
  elseif ~exist(temp_name, 'file') && ~exist(final_name, 'file')
    %xruan 02/21/2016
    % change to slurm job id setting
    % job_id = regexprep(getenv('PBS_JOBID'), '[\r\n]', '');
    job_id = regexprep(getenv('SLURM_JOB_ID'), '[\r\n]', '');
    jid_flag = true;

    if isempty(job_id)
        jid_flag = false;
        p_id = char(num2str(feature('getpid')));
    end
    fid = fopen(temp_name, 'w');
    if jid_flag
        fwrite(fid, sprintf('jid\n%s', job_id));
    else
        fwrite(fid, sprintf('pid\n%s', p_id));
    end
    fclose(fid);    
    can_start = true;
    final_exists = false;
    temp_exists = true;
  elseif exist(temp_name, 'file') && ~exist(final_name, 'file')
    can_start = false;
    final_exists = false;
    temp_exists = true;
  else
    can_start = false;
    final_exists = true;
    temp_exists = true;
  end
  

% % Create unique name using job ID and hostname and use the method described in <http://www.dwheeler.com/secure-programs/Secure-Programs-HOWTO/avoid-race.html>:
% hostname = regexprep(getenv('HOSTNAME'), '[\r\n]', '');
% if isempty(hostname)
% hostname = char(getHostName(java.net.InetAddress.getLocalHost()));
% end

end

