function [z,zmean,zsdev] = tz_zscore(d,dmean,dsdev)
%TZ_ZSCORE Obsolete.
%
%See also ML_ZSCORE

%function [z,zmean,zsdev] = tz_zscore(d,dmean,dsdev)
%OVERVIEW
%   z score transformation
%PARAMETERS
%   d - input data
%   dmean - specified mean
%   dsdev - specified standard deviation
%RETURN
%   z - normalized data
%   zmean - mean of the input data
%   zsdev - std of the input data
%DESCRIPTION
%   
%HISTORY
%   16-May-2005 Initial write TINGZ
%SEE ALSO
%   

error(tz_genmsg('of','tz_zscore','ml_zscore'))

n = size(d,1) ;

if ~exist('dmean','var')
    dmean=[];
end
if isempty(dmean);
    zmean=mean(d,1);
else
    zmean=dmean;
end

if ~exist('dsdev','var')
    dsdev=[];
end
if isempty(dsdev)
    zsdev=std(d,0,1);
    %check constant index
    constidx=tz_constidx(d);
    if ~isempty(constidx)
        zsdev(constidx)=1;
    end
else
    zsdev=dsdev;
end

z = (d-(ones(n,1)*zmean)) ./ (ones(n,1)*zsdev);
