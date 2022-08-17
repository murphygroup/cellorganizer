function [can_start, final_name, final_exists, temp_name] = chunk_start( fname, final_extension, should_block)
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
  %  Tests:
  %  [can_start, final_name, final_exists] = chunk_start('.', 'atomic_test')
  %  dir('./*atomic_test*')
  %  [can_start, final_name, final_exists] = chunk_start('.', 'atomic_test')
  %  dir('./*atomic_test*')
  %  a = pi; save(final_name, 'a'); chunk_finish('.', 'atomic_test');
  %  dir('./*atomic_test*')
  %  [can_start, final_name, final_exists] = chunk_start('.', 'atomic_test')
  %  dir('./*atomic_test*')
  %
  % 2011-12-30 tebuck: Adding job id, if applicable, for greater
  % differences between cluster jobs without hoping for clock or
  % hostname differences.
  % 2012-11-03 tebuck: added an initial file existence check for quick fail.
  % 2013-04-20 tebuck: making final_extension argument optional, adding fourth to sleep until a computing file becomes available.
  % 2013-05-19 tebuck: modifying process so temporary file creation/locking is (hopefully) atomic, jobs do not get done multiple times, and we do not end up with corrupted final files.
  % 2015-01-28 icaoberg: modified code to use matlab's builtin functions instead of depending on system call
  % 2016-02-22 xruan: update enviromemntal variable for job id from PBS to SLURM
  % 2019-09-28 tebuck: do not search for files in unintended directories

  % Process options
  if nargin < 2 || isempty(final_extension)
      final_extension = '.mat';
  end
  if nargin < 3 || isempty(should_block)
      should_block = false;
  end
  
  [fname_path, fname_name, fname_ext] = fileparts(fname);
  if isempty(fname_path)
      fname = fullfile(pwd(), [fname_name, fname_ext]);
  end
  
  temp_name = [ fname, '.tmp'];
  final_name = [ fname, final_extension];

  if ~exist(temp_name, 'file') && exist(final_name, 'file')
    can_start = false;
    final_exists = true;
    return;
  end
  
  %icaoberg 1/28/2015
  %Check which os we're using and adjust the stat function accordingly 
  if ismac()
      statstr = 'stat -L -f"%l" ';
  elseif isunix()
      statstr = 'stat --format="%h" ';
  end
  
  final_exists = false;
  % Create unique name using job ID and hostname and use the method described in <http://www.dwheeler.com/secure-programs/Secure-Programs-HOWTO/avoid-race.html>:
  hostname = regexprep(getenv('HOSTNAME'), '[\r\n]', '');
  if isempty(hostname)
      try
        hostname = char(getHostName(java.net.InetAddress.getLocalHost()));
      catch
          [~, hostname] = system('hostname');
      end
  end
  job_id = regexprep(getenv('SLURM_JOB_ID'), '[\r\n]', '');
  if isempty(job_id)
    job_id = num2str(feature('getpid'));
  end
  unique_name = [ fname, '_', hostname, '_', job_id, '.tmp'];
  current_command = ['touch "', unique_name, '"'];
  [status, output] = system(current_command);
  stat_command = [statstr '"', unique_name, '"'];
  [status, stat_output] = system(stat_command);
  current_command2 = ['link "', unique_name, '" "', temp_name, '"'];
  [status, output] = system(current_command2);
  [status, stat_output2] = system(stat_command);
  
  work_being_done = false;
  
  % hard_link_count = str2num(stat_output)
  % hard_link_count2 = str2num(stat_output2)
  % When code is pasted into the command line, it is taken up by Matlab as if it were part of the stdout of system, so just read the first integer instead of the above two lines:
  hard_link_count = strread(stat_output, '%d', 1);
  hard_link_count2 = strread(stat_output2, '%d', 1);
  % if isempty(hard_link_count) || isempty(hard_link_count2)
    % keyboard
  % end
  % whos hard_link_count hard_link_count2
  if hard_link_count2 == 2 && hard_link_count == 1
    % Go ahead with work, we successfully created the first .tmp file:
    can_start = true;
  else
    % Someone else, possibly another matlabpool worker on the same machine and job, has already linked this or another file, do not start:
    can_start = false;
    work_being_done = true;
  end
  
  % Unlink the unique file, stat_output2 should be 1 if creating a new link fails:
  current_command = ['unlink "', unique_name, '"'];
  [status, output] = system(current_command);
  
  % final_exists, can_start, pause
  
  while true
    % If another job is doing this work but is unfinished and we should wait until the work finishes, just pause and keep checking again if work and writing are both finished:
    if work_being_done && should_block && exist(temp_name, 'file')
      pause(1);
      continue;
    else
      break;
    end
  end
