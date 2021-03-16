function pts = tz_genshape_ac(shapemodel,ncom)
%TZ_GENSHAPE_AC Generate a shape by active shape model.
%   PTS = TZ_GENSHAPE_AC(SHAPEMODEL,NCOM) returns a shape generated from
%   the active shape model SHAPEMODEL. NCOM is the number of components for
%   synthesis.
%   
%   See also

%   26-May-2005 Initial write  T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

ss = shapemodel.avgshape(:)+shapemodel.pcvec(:,1:ncom)* ...
    mvnrnd(zeros(1,ncom),shapemodel.covmat(1:ncom,1:ncom),1)';
pts=[ss(1:360),ss(361:end)];

tz_showpts_2d(pts,'ln',1);
axis('equal')