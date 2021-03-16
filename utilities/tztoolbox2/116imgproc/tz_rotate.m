function pts2 = tz_rotate(pts,theta)
%TZ_ROTATE
%   PTS2 = TZ_ROTATE(PTS,THETA)
%   
%   See also

%   03-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

cost = cos(theta);
sint = sin(theta);
rotateMatrix = [cost,sint;-sint,cost];

pts2 = rotateMatrix*pts;