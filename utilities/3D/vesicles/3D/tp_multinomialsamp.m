function idx = tp_multinomialsamp(P,n)
% Returns index sampled from an multinomial distribution

% Author: Robert F. Murphy and Devin Sullivan
% Edited: Ivan E. Cao-Berg 
%
% Copyright (C) 2011-2012 Murphy Lab
% Lane Center for Computational Biology
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

%Edited by: Devin Sullivan & Bob Murphy 3/1/2012
%   - fixed problem with round off due to large number of points 

P = P/sum(P);

%F = cumulative distribution of probabilties
F = cumsum(P);

%random numbers need to be adjusted for roundoff errors in the cumulative
%sum because there are so many points in P
r = rand(n,1)*F(end);

idx=zeros(length(r),1);
for i = 1:length(r)
     idx(i) = find(r(i)<F,1)+1;
end
