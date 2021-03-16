function subset_idx=tz_sda(features,combclass)
%TZ_SDA Obsolete. See ML_SDA.
%function subidx=tz_sda(features)
%
%OVERVIEW
%   perform sda
%PARAMETERS
%   features - mcf features
%RETURN:
%   subset_idx - index of selected features
%DESCRIPTION:
%   stepwise discriminant analysis algorithm
%
%HISTORY:
%   21-OCT-2004 Initial write TINGZ

warning(tz_genmsg('of','tz_sda','ml_sda'));

if nargin>1
    features=ml_combfeats2mcf(features,combclass);
end

%subidx=mv_stepdisc( features, '~/tmp/sda')
savepath='~/tmp/sda';
hostname = 'linux.andrew.cmu.edu'; % A linux host
tempdir = '~/tmp';

% Prepare input data for SAS
sasinput = ml_fopentemp(tempdir);
sasinput_fullpath = [tempdir '/' sasinput];
ml_make_sasstepdiscinp( features, sasinput_fullpath);

% Prepare command file for SAS
[f, sascmd] = ml_fopentemp(tempdir);
sascmd_fullpath = [tempdir '/' sascmd];
fprintf(f,'%s\n',...
	'PROC IMPORT OUT= WORK.FEATMATRIX',...
	['  DATAFILE= "' sasinput '"'],...
	'  DBMS=CSV REPLACE;',...
	'  GETNAMES=NO;',...
	'  DATAROW=1; ',...
	'RUN;',...
	'proc stepdisc data=FEATMATRIX;',...
	'  class VAR1;',...
	'run;'	);
fclose(f);

% package the data and command files to a tarball
% Ship the data and command files to an andrew unix machine
% because SAS is installed there (I tried to install on our Linux
% machines but was unsuccessful),Run SAS, and ship the results back.
tarcmd = ['tar -C ' tempdir ' -c --bzip -f - ' sasinput ' ' sascmd];
sshcmd = ['ssh -o StrictHostKeyChecking="no" ' hostname];
untarcmd = ['tar -C /tmp -x --bzip -m -f -'];
runsas = ['cd /tmp; sas ' sascmd];
sasoutput = [sascmd '.lst'];
retrieve = ['cat < ' sasoutput];
removetemp = ['rm -f ' sasinput ' ' sascmd '*'];

cmd = [tarcmd ' | ' ...
       sshcmd ' "' untarcmd '; ' runsas '; ' retrieve '; ' removetemp '" ' ...
       '> ' savepath];
fprintf(1,'Andrew unix password: ');
disp(cmd)
keyboard
[status, output] = unix( cmd);
if( status ~= 0), status, output, error(['Error writing tarball for stepdisc']);end

% Delete temporary files
unix(['rm -f ' sasinput_fullpath]);
unix(['rm -f ' sascmd_fullpath]);

% Read the output file from SAS
%disp(savepath)
subset_idx = ml_read_sasstepdiscout( savepath);
