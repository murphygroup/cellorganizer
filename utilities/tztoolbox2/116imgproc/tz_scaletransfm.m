function m2 = tz_scaletransfm(m,scale)
%TZ_SCALETRANSFM Scale transformation matrix.
%   M2 = TZ_SCALETRANSFM(M,SCALE)
%   
%   See also

%   25-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

m2 = [m(1,1) scale(1)/scale(2)*m(1,2) scale(1)*m(1,3);
    m(2,1)*scale(2)/scale(1) m(2,2) m(2,3)*scale(2);
    m(3,1)/scale(1) m(3,2)/scale(2) m(3,3)];
