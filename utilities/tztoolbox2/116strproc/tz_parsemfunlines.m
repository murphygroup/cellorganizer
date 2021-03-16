function [s,n] = tz_parsemfunlines(funlines,igr)
%TZ_PARSEMFUNLINES Parse matlab function declaration.
%   S = TZ_PARSEMFUNLINES(FUNLINES) returns a structure of function
%   declaration by parsing string cell array FUNLINES, in which each
%   element is taken as a line in function declaration. 
%   S has the following fields:
%       input: a array of input arguments
%       output: a cell array of returned arguments
%       funname: function name
%
%   TZ_PARSEMFUNLINES(FUNLINES,IGNORE) will ignore keyword 'function'
%   in FUNLINES when IGNORE is not 0.
%
%   [S,N] = TZ_PARSEMFUNLINES(...) also returns number of lines N
%   in function declaration.
%
%   Example:
%       [s,n] = tz_parsemfunlines({'[a,b]=test(c,...','d)'},1)
%       returns the following values:
%           s.input={'c';'d'}
%           s.output={'a';'b'}
%           s.funname='test'
%           n=2

%   15-Apr-2005 Initial write TINGZ
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin<1
    error('1 or 2 arguments are required');
end

if nargin<2
    igr=0;
end

if ~iscell(funlines)
    error('The first input must be a cell array');
end

warnmsg=[];
s.input={};
s.output={};
s.funname='';

%extract function definition
n=0;

%In case the first line is 'function...'
if(strcmp(funlines{1},'function...'))
    funlines{1}='function ...';
end

%count lines of function of declaration
for i=1:length(funlines)
    n=n+1;
    tmpline=funlines{i};
    tmpline=strrep(tmpline,' ','');
    if length(tmpline)<3
        break;
    end
    
    if ~strcmp(tmpline(end-2:end),'...')
        break;
    end    
end

%combines all lines of function declaration together
fundef=tz_cell2str(funlines(1:n),'');
fundef=strrep(fundef,'...','');

%splitted into strings by space
fundefds=strread(fundef,'%s','delimiter',' ');


isfun=1;

%searching for keyword 'function'
if isempty(fundefds)
    isfun=0;
else   
    if ~strcmp(fundefds{1},'function')
        if(igr==0)
            isfun=0;  
        end
    else
        fundefds(1)=[];
    end
end

if(isfun==0) | isempty(fundefds)
    s.warnmsg{1}='not a function';
    return;
end


fundef=tz_cell2str(fundefds,'');

%searching for keyword '='
fundefds=strread(fundef,'%s','delimiter','=');
if length(fundefds)==1
    s.output={};
else
    output=fundefds{1};
    output=tz_strdel(output,'[]');
    s.output=strread(output,'%s','delimiter',',');
    fundef=tz_cell2str(fundefds(2:end),'');
end


%extracting function name
fundefds=strread(fundef,'%s','delimiter','(');
if length(fundefds)==1
    s.funname=fundefds{1};
    s.input={};
    return
else
    s.funname=fundefds{1};
    input=fundefds{2};
    input=tz_strdel(input,')');
    s.input=strread(input,'%s','delimiter',',');
end