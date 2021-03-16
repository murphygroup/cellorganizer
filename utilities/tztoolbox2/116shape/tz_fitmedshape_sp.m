function spmodel = tz_fitmedshape_sp(axln,dist,knots,k)
%TZ_FITMEDSHAPE_SP Fit a median axis shape model by B-spline (Obsolete).
%   SPMODEL = TZ_FITMEDSHAPE_SP(AXLN,DIST,KNOTS,K) returns the spline
%   medial axis desription of the shape with medial axis AXLN and width
%   DIST. KNOTS and K are number of knots and order for the splines. 
%   
%   Notice: The returned structure is not 'mxp' structure.
%   
%   See also

%   26-May-2005 Initial write  T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 4
    error('Exactly 4 arguments are required')
end

y=dist;
len=length(y);
spmodel.len=len;
x=(0:(len-1))/(len-1); 

%figure
sp=spap2(knots,k,x,y);
    
if knots>1
    spmodel.dspara=[sp.knots(k+1:k+knots-1),sp.coefs];
else
    spmodel.dspara=[sp.coefs];
end
spmodel.dssp=sp;

fnplt(sp);
hold on
plot(x,y,'rx');
hold off
drawnow

y=axln(:,2)';
sp=spap2(knots,k,x,y)
if knots>1
    spmodel.lnpara=[sp.knots(k+1:k+knots-1),sp.coefs];
else
    spmodel.lnpara=[sp.coefs];
end
spmodel.lnsp=sp;

figure
fnplt(sp);

hold on
plot(x,y,'rx');
hold off
drawnow

