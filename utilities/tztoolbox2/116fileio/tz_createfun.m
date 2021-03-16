function tz_createfun(dirname,filename,inputs,outputs,optionalinputs, ...
    openmode)
%TZ_CREATEFUN Automatic generate a template code for matlab function.
%   TZ_CREATEFUN(DIRNAME,FILENAME) generates the file for a fuction
%   with the name string FILENAME under directory DIRNAME. The function
%   has no input and output. The comments have 'tz' style. If FILENAME
%   is a cell array, FILENAME will be parsed as matlab function
%   definition.
%   
%   TZ_CREATEFUN(DIRNAME,FILENAME,INPUTS) is for a function with
%   input arguments INPUTS, which is a cell array of names of arguments
%   necessary for the function. The comments have 'tz' style.
%
%   TZ_CREATEFUN(DIRNAME,FILENAME,INPUTS,OUTPUTS) is for function
%   with necessary input argument INPUTS and output argument
%   OUTPUTS. The comments have 'tz' style.
%   
%   TZ_CREATEFUN(DIRNAME,FILENAME,INPUTS,OUTPUTS,OPTIONALINPUTS)
%   also consider optional input arguments for the function.
%   The comments have 'mt' style.
%
%   TZ_CREATEFUN(DIRNAME,FILENAME,INPUTS,OUTPUTS,OPTIONALINPUTS,OPENMODE)
%   speicifies how to open the created file by OPENMODE:
%       0 : do not open
%       1 : open by matlab editor. This is the default value.
%       2 : open by emacs
%
%   Example:
%       tz_createfun('./','test',{'i1','i2'},{'o1'},{'i3'})
%       creates function o1=test(i1,i2,i3) under current
%       directory. The editor will automatically open the created
%       file, which is test.m in this case.
%       tz_createfun('./','test',{'i1','i2'},{'o1'}) is the same
%       as tz_createfun('./',{'o1=test(i1,i2)'}).
%
%   See also TZ_GENCOMMENTS TZ_GENARGCHECK

%   12-Aug-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('At least 2 arguments are required')
end

if nargin<5
    style='tz';
    optionalinputs={};
else
    style='mt';
end

narg = nargin;

if iscell(filename)
    [s,n]=tz_parsemfunlines(filename,1);
    filename=s.funname;
    inputs=s.input;
    outputs=s.output;
    narg = 4;
end

filepath=[dirname '/' filename '.m'];

isexist = 0;

if ~isempty(which(filename))
    warning(['Name confilct: Another function has the same name: ' which(filename)]);
end

if ~exist(filepath,'file')   
    if narg < 4
        outputs={};
        optionalinputs={};
    end
    
    if narg < 3
        inputs={};
        outputs={};
        optionalinputs={};
    end
    
    fid=fopen(filepath,'w');
    
    funtitle=['function ' tz_funcell2str(filename, ...
        {inputs{:},optionalinputs{:}},outputs)];
    
    if strcmp(tz_getuser,'tingz')
        author = 'T. Zhao';
    else
        author = tz_getuser;
    end
    
    comments=tz_gencomments(filename,inputs,outputs, ...
        optionalinputs,style,author);
    
    fprintf(fid,'%s\n',funtitle);
    fprintf(fid,'%s\n',comments);
    fprintf(fid,'\n');
    fprintf(fid,'%s',tz_genargcheck(length(inputs),length(optionalinputs))); 
    
    fclose(fid);
else 
    warning('function alread exists');
    isexist = 1;
end

if ~exist('openmode','var')
    openmode = 1;
end

switch openmode
    case 0
        if ~isexist
            disp(['function ' filename ' created']);
        end   
    case 1
        open(filename);
    case 2
        tz_open(filename);
    otherwise
        error('Unrecognized open mode');
end


