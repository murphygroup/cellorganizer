function [hd, lhd, rhd, d_1, d_2] = hausdorff_distance_by_nearest_neighbor_search(X_s, X_t)
% compute hausdorf distance
% 04/10/2018 fix bug for N_ptn
% 09/20/2018 use nearest neighbor search to calculate hausdorff distance

%
% Author: Xiongtao Ruan
%
% Copyright (C) 2012-2018 Murphy Lab
% Computational Biology Department
% School of Computer Science
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



[~, d_1] = knnsearch(X_t, X_s);
lhd = max(d_1);
[~, d_2] = knnsearch(X_s, X_t);
rhd = max(d_2);

hd = max(lhd, rhd);

end