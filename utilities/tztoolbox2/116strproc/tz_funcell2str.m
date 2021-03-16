function funtitle = tz_funcell2str(funname,inputs,outputs)
%TZ_FUNCELL2STR Build a full function name from parts
%   TZ_FUNCELL2STR(FUNNAME) returns FUNNAME.
%   
%   TZ_FUNCELL2STR(FUNNAME,INPUTS) returns a function name with input
%   arguments stored in the cell array INPUTS.
%
%   TZ_FUNCELL2STR(FUNNAME,INPUTS,OUTPTS) also adds the output arguments
%   OUTPUTS into the function name. Both INPUTS and OUTPUTS 
%   are cell arrays.
%
%   Example:
%       tz_funcell2str('testfun',{'arg1','arg2'}) returns
%       'testfun(arg1,arg2)'
%
%       tz_funcell2str('testfun',{'arg1','arg2'},{'o1','o2','o3'}) 
%       returns '[o1,o2,o3] = testfun(arg1,arg2)'    

%   03-Aug-2005 Initial write TINGZ
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin<1
    error('At least one argument is required');
end

if ~exist('outputs','var')
    outputs={};
end

if ~exist('inputs','var')
    inputs={};
end

if ~iscell(inputs)
    error('The second argument must be a cell array');
end

if ~iscell(outputs)
    error('The third argument must be a cell array');
end

if ~isempty(outputs)
    if(length(outputs{1}) > 0)
        if length(outputs)==1
            outform=[outputs{1}];
        else
            outform=['[' tz_cell2str(outputs) ']'];
        end
        outform=[outform ' = '];
    else
        outform=[];
    end
else
    outform=[];
end

if ~isempty(inputs)
    if length(inputs{1})>0
        inform=['(' tz_cell2str(inputs) ')'];
    else
        inform=[];
    end
else
    inform=[];
end

funtitle=[outform funname inform];