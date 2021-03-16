function m = ml_wmoment(x,ws,order)

%ML_WMOMENT calculates weighted central moments
%   M=ML_WMOMENT(X,WS,ORDER) the central moment of X specified by the
%   ORDER. WS is the weights. Both X and WS must be a vector and have the
%   same size. However, The first order will return moment, the mean,
%   rather than the central moment.

% 24-Apr-2005 Initial write TINGZ
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


if order==1
    m=sum(x.*ws)/sum(ws);
    return
end

mu=ml_wmoment(x,ws,1);
m=sum((x-mu).^order.*ws)/sum(ws);
