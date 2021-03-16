function [compartment_stats] = compartment_percell(segdna, segcell, improt, param)
%This code computes the per-cell compartmental statistics
%Compartments currently used are - plasma membrane('mem'), 
%cytoplasm('cyto'),and nuclear('nuc')

% Author: Devin Sullivan 6/11/13
%
% Copyright (C) 2007-2013 Murphy Lab
% Carnegie Mellon University
%
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
% For additional information visit http://murphylab.web.cmu.edu/ or
% send email to murphy@cmu.edu

%get the current seg results 


% segresults = [param.preprocessingFolder filesep 'cell.mat'];
% load(segresults,'segdna','segcell','top_slice','bot_slice');
% param.topslice = top_slice;
% param.botslice = bot_slice;
%Define the compartments 

[indexedimg,compartmentlist] = DefineCompartments(segcell,segdna,param);

%Get the fluorescence in each compartment
if ismember('prot', compartmentlist)
    [totfluor,compfluor,propfluor] = ...
        CompartmentProportions(improt,indexedimg,compartmentlist,param);

        compartment_stats.totfluor = totfluor;
    compartment_stats.compfluor = compfluor;
    compartment_stats.propfluor = propfluor;
end

compartment_stats.indexedimg = indexedimg;
compartment_stats.compartmentlist = compartmentlist;


