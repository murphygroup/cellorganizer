function [pts2,coef1,coef2]=tz_ftshape(pts,coefnum,isshow)
%TZ_FTSHAPE Shape description by Fourier transform.
%   PTS2 = TZ_FTSHAPE(PTS,COEFNUM,ISSHOW) returns an array of points that
%   are reconstruction from the first COEFNUM components of Fourier
%   transform. If ISSHOW is 1, the resconstructed points will be plotted.
%   
%   [PTS2,COEF1,COEF2] = TZ_FTSHAPE(...) aksi returns the coefficients of
%   Fourier transform.
%   
%   See also TZ_WAVSHAPE

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function [pts2,coef1,coef2]=tz_ftshape(pts,coefnum)

ft1=fft(pts(:,1));
ft2=fft(pts(:,2));

ft1(coefnum+1:end)=0;
ft2(coefnum+1:end)=0;

coef1=ft1;
coef2=ft2;

pts2=[abs(ifft(coef1)),abs(ifft(coef2))];

if isshow
    plot(pts2(:,1),pts2(:,2));
    axis('equal')
end