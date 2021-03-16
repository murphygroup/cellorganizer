function [testnames,missed,assigned]=ml_missclassnames(names,testidx,confmatout)
% ML_MISSCLASSNAMES - Returns the file names of misclassified samples.
%
% [MISSED TRUTH]=ML_MISSCLASSNAMES(NAMES,TESTIDX,CONFMATOUT)
%       Outputs:
%
%       Inputs:
%
%
%  M. Boland - 01 May 1999
%

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

% $Id: ml_missclassnames.m,v 1.2 2006/06/27 13:33:47 tingz Exp $

if (~iscell(names) | ~iscell(testidx) | ~iscell(confmatout))
  error('NETOUT, CLASS, and CONFMATOUT must be cell arrays') ;
end

if (length(names)<1)
  error('NAMES contains no elements') ;
end

numc = length(names) ;
numtrials = length(confmatout) ;

testnames={} ;
missed={} ;
assigned=[] ;

for i=1:numtrials 
  thesetestnames={} ;
  for j=1:numc
    thesetestnames = [thesetestnames ; names{j}(testidx{j}{i})] ;
  end
  testnames = [testnames ; thesetestnames] ;
  missed = [missed ; thesetestnames(confmatout{i}.missed.index)] ;
  assigned = [assigned ; confmatout{i}.missed.assigned'] ;
end

testnames = sort(testnames) ;
[missed, sortidx] = sort(missed) ;
assigned = assigned(sortidx) ;

