function tz_updaterecord(events,etype,time)

%function tz_updaterecord(events)
%
%OVERVIEW:
%   Record history to the file update.txt
%PARAMETERS:
%   events - a string containing events for saving
%   etype - the type of event
%   time - the time at which the events happen
%          default is current time
%RETURN:
%
%DESCRIPTION
%   This function is for recording events such as saving results, modifying files, etc
%   The etype could be 'update','results','chname','other' etc

if ~exist('time','var')
    time=num2str(round(clock));
end

load local.mat

fid=fopen(['~/matlab/' machine '/update.txt'],'a');

updatetime=round(clock);

fprintf(fid,['@ ' num2str(round(updatetime)) '\n']);
if ~isempty(etype)
    fprintf(fid,['!' etype '\n']);
end
fprintf(fid,[events '\n']);

fclose(fid);




