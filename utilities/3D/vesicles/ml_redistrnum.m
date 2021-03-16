function n=ml_redistrnum(N,k)
%ML_REDISTRNUM Evenly separate a number into several numbers.
%   N = ML_REDISTRNUM(M,K) separated the number M into a vector with the
%   length K evenly. The sum of N is M.
%   
%   Example:
%       ml_redistrnum(10,4) returns [3 3 2 2].

%   ??-???-???? Initial write T. Zhao
%   05-NOV-2004 Modified T. Zhao
%       - add comments
%   22-MAR-2005 Modified T. Zhao
%       - debug k==1
%   Copyright (c) Murphy Lab, Carnegie Mellon University

% Copyright (C) 2007  Murphy Lab
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

if nargin < 2
    error('Exactly 2 arguments are required')
end

if k==1
    n=N;
    return
end

avgn=floor(N/k);
for i=1:k
    n(i)=avgn;
end

remain=N-sum(n);
k=1;
while remain>0
    n(k)=n(k)+1;
    remain=remain-1;
    k=k+1;
end

% 
% n(k)=N-sum(n(1:k-1));
% 
% [maxn,maxpos]=max(n);
% [minn,minpos]=min(n);
% 
% while(maxn-minn>1)
%     n(maxpos)=n(maxpos)-1;
%     n(minpos)=n(minpos)+1;
%     [maxn,maxpos]=max(n);
%     [minn,minpos]=min(n);
% end