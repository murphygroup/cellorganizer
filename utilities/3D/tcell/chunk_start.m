function [can_start, final_name, final_exists, temp_exists] = chunk_start(path, fname, final_extension, should_block)
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
  % 2014-05-27 tebuck: shouldn't this be set up so final_name is also temporary and is copied to the final destination only by chunk_finish?
  
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
  end
  
  final_exists = false;
  % Create unique name using job ID and hostname and use the method described in <http://www.dwheeler.com/secure-programs/Secure-Programs-HOWTO/avoid-race.html>:
  hostname = regexprep(getenv('HOSTNAME'), '[\r\n]', '');
  if isempty(hostname)
    hostname = char(getHostName(java.net.InetAddress.getLocalHost()));
  end
  %xruan 02/21/2016
  % change to slurm job id setting
  % job_id = regexprep(getenv('PBS_JOBID'), '[\r\n]', '');
  job_id = regexprep(getenv('SLURM_JOB_ID'), '[\r\n]', '');
  
  if isempty(job_id)
    job_id = char(num2str(feature('getpid')));
  end
  unique_name = [path, '/', fname, '_', hostname, '_', job_id, '.tmp'];
  current_command = ['touch "', unique_name, '"'];
  [status, output] = system(current_command);
  stat_command = ['stat --format="%h" "', unique_name, '"'];
  [status, stat_output] = system(stat_command);
  current_command = ['link "', unique_name, '" "', temp_name, '"'];
  [status, output] = system(current_command);
  [status, stat_output2] = system(stat_command);
  
  work_being_done = false;
  if length(stat_output) > 2 || length(stat_output2) > 2 || length(stat_output) == 0 || length(stat_output2) == 0
      pause(1);
      if length(stat_output) > 2 || length(stat_output) == 0
          if ~exist(unique_name, 'file')
                fclose(fopen(unique_name, 'w'));
          end
          stat_output = '1';
      end

      if length(stat_output2) > 2 || length(stat_output2) == 0
          for i = 1 : 5
              current_command2 = ['link "', unique_name, '" "', temp_name, '"'];
              [status, output] = system(current_command2);
              [status, stat_output2] = system(stat_command);
              if length(stat_output2) == 2 && strcmp(stat_output2(1), '2')
                  break;
              else
                  stat_output2 = '0';
              end
          end
      end

  end
  % hard_link_count = str2num(stat_output)
  % hard_link_count2 = str2num(stat_output2)
  % When code is pasted into the command line, it is taken up by Matlab as if it were part of the stdout of system, so just read the first integer instead of the above two lines:
  hard_link_count = strread(stat_output, '%d', 1);
  hard_link_count2 = strread(stat_output2, '%d', 1);
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
    temp_exists = exist(temp_name, 'file') > 0;
    if work_being_done && should_block && temp_exists
      pause(1);
      continue;
    else
      final_exists = exist(final_name, 'file') > 0;
      break;
    end
  end

