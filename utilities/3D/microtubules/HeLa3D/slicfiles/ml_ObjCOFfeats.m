function [meandist, stddist, maxmindist] = ml_ObjCOFfeats( dists)

% [MEANDIST, STDDIST, MAXMINDIST] = ML_OBJCOFFEATS( DISTS)
%
% Calculates the three features relating objects COF-s to a
% central COF:
% 1) The average object distance to central COF
% 2) The stddev of object distances to central COF
% 3) The Max/Min ratio of object distances to central COF
% DISTS is a vector of object distances to COF

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

meandist = mean( dists); % Average object dist to COF
stddist = std( dists); % Std.dev. of object dist to COF
if( min( dists) > 0)
    %Ratio of smallest to largest obj to COFdist
    maxmindist = max( dists) / min( dists); 
else
    maxmindist = 0;
end
