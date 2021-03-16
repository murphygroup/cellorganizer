function tz_parseupdate(option,varargin)

%TZ_PARSEUPDATE Parse updating files.

%function tz_parseupdate(option)
%   
%OVERVIEW:
%   parse the updating file
%PARAMETERS:
%   option - the behavior of the function
%   varargin - other arguments
%       - arbitary number of arguments
%RETURN:
%
%DESCRIPTION:
%
%HISTORY:
%   ??-???-2004 Initial write TINGZ
%   02-NOV-2004 Modified TINGZ
%       - add comments

load local.mat

fid=fopen(['~/matlab/' machine '/update.txt'],'r');
time=[];
etype=[];
numtime=[];
isevent=0;

while 1
    switch(option)
    case 'disp'
        isdisp=1;
        for i=1:length(varargin)
            if strcmp(varargin{i},'type')
                if length(etype)==length(varargin{i+1})
                    if ~all(etype==varargin{i+1})
                        isdisp=[isdisp 0];
                    end
                end
            end
            
            if strcmp(varargin{i},'time')
                if ~isempty(numtime)
                    if isempty(varargin{i+2})
                        varargin{i+2}=Inf;
                    end
                    if ~(tz_timecmp(numtime,varargin{i+1})>=-1 & tz_timecmp(numtime,varargin{i+2})<=1)
                        isdisp=[isdisp 0];
                    end
                end
            end
        end
        
        if all(isdisp==1) & isevent==1
            disp(['@' time]);
            disp(event);
        end
        
    end
    
    tline = fgetl(fid);
    
    if ~ischar(tline), break, end
    
    isevent=0;
    switch tline(1)
    case '@'
        time=tline(2:end);
        numtime=str2num(time);
    case '!'
        etype=tline(2:end);
    otherwise
        event=tline;
        isevent=1;
    end
end

fclose(fid);