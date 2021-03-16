function peaks=tz_findpeak(x,option)
%TZ_FINDPEAK Find peaks of a curve.
%   PEAKS = TZ_FINDPEAK(X,OPTION) returns local maxima or minima or both
%   of the curve specified by 1D array X. PEAKS has two column. Its first
%   column contains indices and second column contains corresponding
%   values. OPTION can be specified as follows:
%       'pk' - local maxima
%       'vl' - local minima
%       'pv' - both
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

dx=x-[x(end) x(1:end-1)];
dx=dx>0;
px=dx-[dx(2:end) dx(1)];

switch option
case 'pk'
    peaks(1,:)=find(px>0); 
case 'vl'
    peaks(1,:)=find(px<0);
case 'pv'
    peaks(1,:)=find(px~=0);
end

peaks(2,:)=x(peaks(1,:));
peaks=peaks';
