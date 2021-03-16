function [summary]=ml_mlpclasssummary(confmats)
% ML_MLPCLASSSUMMARY - Summarize the output from MB_MLPCONFMAT
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

% $Id: ml_mlpclasssummary.m,v 1.2 2006/06/27 13:33:47 tingz Exp $

if(~iscell(confmats))
  error('CONFMATS must be a cell array in which each element contains the output from ml_mlpconfmat') ;
end

if(length(confmats)<1)
  error('CONFMATS has no elements') ;
end

cmat = zeros(size(confmats{1}.cmat)) ;
crates = [] ;

for i=1:length(confmats)
  crates = [crates ; confmats{i}.crate(:,1:2)] ; 
  cmat = cmat+confmats{i}.cmat ;
end

summary.confusion = cmat ./ (sum(cmat')' * ones(1,size(cmat,2))) ;

cmat_nounk = cmat(:,1:(end-1)) ;
summary.confusion_nounk = cmat_nounk ./ (sum(cmat_nounk')' * ...
                   ones(1,size(cmat_nounk,2))) ;

summary.Pc_mean = [mean(diag(summary.confusion)) ...
                   mean(diag(summary.confusion_nounk))] ;

summary.Pc_var = var(crates) ;

