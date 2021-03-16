function [pts2,coef1,coef2]=tz_wavshape(pts,level,isshow)
%TZ_WAVSHAPE Wavelet shape description.
%   PTS2 = TZ_WAVSHAPE(PTS,LEVEL,ISSHOW) returns an array of points that
%   are reconstruction from the approximate coefficients of  wavelet 
%   trasform at LEVEL level. If ISSHOW is 1, the resconstructed points 
%   will be plotted.
%   
%   [PTS2,COEF1,COEF2] = TZ_WAVSHAPE(...) also returns the wavelet
%   transform coefficients.
%   
%   See also TZ_FTSHAPE

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('Exactly 3 arguments are required')
end

if level==0
    pts2=pts;
    coef1=pts(:,1);
    coef2=pts(:,2);
else
    
    [c1,l1]=wavedec(pts(:,1),level,'db4');
    [c2,l2]=wavedec(pts(:,2),level,'db4');
    
    coef1=c1(1:l1(1));
    coef2=c2(1:l2(2));
    
    pts2=[waverec(coef1,l1(1:2),'db4'),waverec(coef2,l2(1:2),'db4')];
end
if isshow
    plot(pts2(:,1),pts2(:,2));
    axis('equal')
end