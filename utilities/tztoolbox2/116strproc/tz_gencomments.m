function comments = ...
    tz_gencomments(functionname,inputs,outputs,optionalinputs,style,author)
%TZ_GENCOMMENTS Generates matlab comments automatically.
%   TZ_GENCOMMENTS(FUNCTIONNAME,INPUTS,OUTPUTS,OPTIONALINPUTS,STYLE)
%   returns a string of comment template of a function with the name
%   FUNCTIONNAME. INPUTS is a string array of required input arguments. It
%   is empty if the function can take no argument. OUTPUTS is a string
%   array of return values. It is empty if there is no output argument.
%   OPTINALINPUTS is a string array of optional input arguments. It is
%   empty if all input arguments are necessary. STYLE specifies the style
%   of commenting. There are two choices, 'mt' and 'tz'. 'mt' is
%   recommended because it is has standard matlab format.
%   INPUTS,OUTPUTS or OPTIONALINPUTS can also be a string, in which 
%   variables are separated by ',' or space.
%
%   TZ_GENCOMMENTS(FUNCTIONNAME,INPUTS,OUTPUTS,OPTIONALINPUTS,STYLE,AUTHOR)
%   specifies the author of the script. Default is 'T. Zhao'.
%
%   Example:
%   for fucntion [pvalue,ts] = ml_ht2test2(x,y,pooled), we can use
%     tz_gencomments('ml_ht2test2','x,y','pvalue,ts','pooled','mt') 
%   or
%     tz_gencomments('ml_ht2test2',{'x','y'},{'pvalue','ts'},'pooled','mt') 

%   03-Aug-2005 Initial write T. Zhao
%   18-Sep-2005 Modified T. Zhao
%       - support strings for inputs,outputs and optionalinputs
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin<5
    error('5 or 6 arguments are required');
end
   
if ~exist('author','var')
    author = 'T. Zhao';
end

varSeparators = ', ';
if isstr(inputs)
    inputs = tz_strtok(inputs,varSeparators);
end

if isstr(outputs)
    outputs = tz_strtok(outputs,varSeparators);
end

if isstr(optionalinputs)
    optionalinputs = tz_strtok(optionalinputs,varSeparators);
end

if ~iscell(inputs) | ~iscell(outputs) | ~iscell(optionalinputs)
    error('INPUTS, OUTPUTS and OPTIONALINPOUTS should all be cell arrays or strings.');
end

if isempty(outputs)
    outputs = {''};
end

indent = '   ';
commentChar = '%';
copyright = 'Murphy Lab';
copyright2 = 'Carnegie Mellon University';

switch(style)
case 'mt'
    %First line of the comments.
    comments = [commentChar funusage(functionname,{},{},style)];
    
    %Comments for no optional arguments.
    comments = addcmtline(comments,[commentChar indent ...
            funusage(functionname,inputs,outputs(1), style)]);

    comments = addcmtline(comments,[commentChar,indent]);
    
    %Comments for optional arguments.
    for i = 1:length(optionalinputs)
        comments = addcmtline(comments,[commentChar indent ... 
                funusage(functionname, { inputs{:},optionalinputs{1:i} }, ...
                outputs(1), style)]);
        comments = addcmtline(comments,[commentChar,indent]);
    end
    
    if(length(outputs) > 1)
        inputs={'...'};
        comments = addcmtline(comments,[commentChar indent ...
                funusage(functionname,inputs,outputs, style)]);
        comments = addcmtline(comments,[commentChar,indent]);
    end
    
    comments = addcmtline(comments,[commentChar indent 'See also']);
    
case 'tz'
    comments = '';
    funtitle=funusage(functionname, { inputs{:},optionalinputs{:} }, ...
                outputs, style);
    comments = addcmtline(comments,[commentChar funtitle]);
    
    keywords=tz_getcmtkeywords;
    comments = addcmtline(comments,[commentChar keywords{1}]);  %overview
    comments = addcmtline(comments,[commentChar indent]);
    comments = addcmtline(comments,[commentChar keywords{2}]);  %arguments
    
    if ~isempty(inputs)
        for i=1:length(inputs)
            comments=addcmtline(comments,[commentChar indent inputs{i} ' - ']);
        end
    else
        comments=addcmtline(comments,[commentChar]);
    end
    
    comments = addcmtline(comments,[commentChar keywords{3}]); %output
    if ~isempty(outputs)
        for i=1:length(outputs)
            comments=addcmtline(comments,[commentChar indent outputs{i} ' - ']);
        end
    else
        comments=addcmtline(comments,[commentChar]);
    end
    
    comments = addcmtline(comments,[commentChar keywords{4}]); %description
    comments = addcmtline(comments,[commentChar indent]);
    
    comments = addcmtline(comments,[commentChar keywords{5} ' ']); %see also
    comments = addcmtline(comments,[commentChar indent]);
otherwise
    error([style ': Invalid style.'])
end

comments = addcmtline(comments,[]); %end of comments for help

if strcmp(style,'tz')
    comments = addcmtline(comments,[commentChar keywords{6}]); %history
end

comments = addcmtline(comments,[commentChar indent date ' Initial write ' ...
                    author]);
year = tz_strtok(date,'-');
year = year{end};
comments = addcmtline(comments,[commentChar indent 'Copyright (c) ' year ...
                    ' ' copyright]);
comments = addcmtline(comments,[commentChar indent copyright2]);
comments = addcmtline(comments,commentChar);
licenseStatements = help('ml_gpl');
licenseStatements = [commentChar ' ' strrep(licenseStatements, ...
               [char(10)],[char(10) commentChar ' '])];
licenseStatements = licenseStatements(1:end-2);
comments = addcmtline(comments,licenseStatements);

%%%%%% subfunction %%%%%%%%%%%%%%%%
function comments = funusage(functionname,inputs,outputs,style)

nRequiredArguments=4;

if nargin~=nRequiredArguments
    error(['Exactly ' num2str(nRequiredArguments) ' arguments are required']);
end

comments=[tz_funcell2str(functionname,inputs,outputs)];

switch(style)
case 'mt'
    comments=[upper(comments)];    
case 'tz'
    comments=['function ' comments];
otherwise
    error('Invalid style.')
end

function comments=addcmtline(comments,cmtline)

comments=tz_addstrline(comments,cmtline);