function Feat=ml_wavefeatures(image);

%FEAT=ML_WAVEFEATURES(IMAGE); constructs a matrix composed of the
%wavelet features for the given image based on detail coefficients from
% 'db4' decomposed 10 levels deep.
%
%It's input is a normalized image

% Copyright (C) 2006  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

Feat=zeros(0,1);

[C,S] = wavedec2(image,10,'db4');

for i = 0 : 9
    chd = detcoef2('h',C,S,(10-i));
    cvd = detcoef2('v',C,S,(10-i));
    cdd = detcoef2('d',C,S,(10-i));
    
    hfeat = sqrt(sum(sum(chd.^2)));
    vfeat = sqrt(sum(sum(cvd.^2)));
    dfeat = sqrt(sum(sum(cdd.^2)));
    
    Feat = [Feat hfeat vfeat dfeat];
end;
