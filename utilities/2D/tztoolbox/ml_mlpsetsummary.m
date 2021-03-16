function [confusion, confusion_nounk, Pc]=ml_mlpsetsummary(setinfo)
% ML_MLPSETSUMMARY - Summarize the output from ML_MLPCONFMAT
%
% M. Boland - 14 Apr 1999

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

% $Id: ml_mlpsetsummary.m,v 1.2 2006/06/27 13:33:47 tingz Exp $

if(~iscell(setinfo))
  error('SETINFO must be a cell array in which each element contains the output from ml_mlptrainstoptestsets') ;
end

if(length(setinfo)<1)
  error('SETINFO has no elements') ;
end

cmat = zeros(size(setinfo{1})) ;

for i=1:length(setinfo)
  cmat = cmat+setinfo{i} ;
end

confusion = cmat ./ (sum(cmat')' * ones(1,size(cmat,2))) ;

cmat_nounk = cmat(:,1:(end-1)) ;
confusion_nounk = cmat_nounk ./ (sum(cmat_nounk')' * ...
                   ones(1,size(cmat_nounk,2))) ;

Pc = [mean(diag(confusion)) mean(diag(confusion_nounk))] ;



