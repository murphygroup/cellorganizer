function t=tz_timecmp(t1,t2)

%function t=tz_timecmp(t1,t2)
%
%OVERVIEW
%   compare time difference
%PARAMETERS:
%   t1 - time 1 num array [year month day hour min sec]
%   t2 - time 2 num array [year month day hour min sec]
%RETURN:
%   t - which one is later
%       2 - t1 is later
%       -2 - t2 is later
%       0 - the same
%       1 - undecided, t1 has higher accuracy
%       1 - undecided, t2 has higher accuracy
%DESCRIPTION:
%   
%HISTORY:
%   ??-???-2004 Initial write
%   03-NOV-2004 Modified TINGZ
%

%?
t1=num2cell(t1);
t2=num2cell(t2);

minacc=min([length(t1),length(t2)]);

t=0;

for i=1:minacc
    if t1{i}>t2{i}
        t=2;
        break
    end
    
    if t1{i}<t2{i}
        t=-2;
        break
    end
end

if t==0
    if length(t1)>length(t2)
        t=1;
    end
    
    if length(t1)<length(t2)
        t=-1;
    end
end

