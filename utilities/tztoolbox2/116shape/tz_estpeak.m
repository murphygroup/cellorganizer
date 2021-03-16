function peaks=tz_estpeak(x,option)
%TZ_ESTPEAK Estimate peaks of a curve.
%   PEAKS = TZ_ESTPEAK(X,OPTION) returns the local maxima or minima or both
%   of the curve estimated from the vector X. See TZ_FINDPEAK for details
%   about OPTION.
%   
%   See also TZ_FINDPEAK

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

y=[x,x];

sp = spap2(45,5,1:length(y),y);

yy=spval(sp,1:length(y));

peak=tz_findpeak(yy,option);

peaks=[peak(2,:);peak(3,:)];

for i=4:size(peak,1)
    
    if peaks(end,1)-peaks(1,1)>length(x)-5
        break;
    end
    peaks=[peaks;peak(i,:)];
end

peaks(end,:)=[];

if any(peaks(:,1)>360)
    peaks(peaks(:,1)>360)=peaks(peaks(:,1)>360)-360;
end

[peaks(:,1),I]=sort(peaks(:,1));
peaks(:,2)=peaks(I,2);
