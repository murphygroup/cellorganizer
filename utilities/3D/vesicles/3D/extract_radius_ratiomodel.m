function [rad_ratio,cellnucratios,nucbottomslices] = extract_radius_ratiomodel(cellcodepath)
% From the cellcodes, extract the radius ratio of nuclear radius over cell
% radius under cylindrical systerm

% Author: Devin Sullivan - adapted from Tao Peng's extract_radius_ratio.m
% Edited: Ivan E. Cao-Berg
%
%%
% June 10 2013 D. Sullivan extract_radius_ratiomodel builds the whole
%                          pattern model using per-cell parameters 
%                          precomputed using extract_radius_ratio.m
% July 9, 2013 G. Johnson  Bug check to make sure rad_ratio lengths are the
%                           same
% July 23, 2013 G. Johnson Make sure no NaN or Inf values are added to the
%                          ratio model
%%
%
% Copyright (C) 2011 Murphy Lab
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

% Constants
%protype = {'DNA','LAM','Nuc','Gia','Gpp','Mit','TfR'};
n = 360;
H = -pi:(2*pi/360):pi;


% Calculate radius ratios and combining them
celllist = ml_dir([cellcodepath filesep 'cellcodes_*']);
load([cellcodepath filesep celllist{1}]);
distratio = zeros(length(celllist),length(rad_ratio)); 
cellnucratios = zeros(length(celllist),1); 
nucbottomslices = zeros(length(celllist),1);

%D. Sullivan 6/10/13 
%Loop through the per-cell results and compile the radius ratios to build a
%model
k = 0;
for i = 1:length(celllist)
    fname = [cellcodepath filesep celllist{i}];
    if exist(fname,'file')
        
        load(fname)
        
        %G. Johnson 7/23/23 added checks to not allow inf or nans
        if length(rad_ratio) == size(distratio,2) && ~any(rad_ratio == inf) && ~any(isnan(rad_ratio))
            k = k + 1;
            distratio(k,:) = rad_ratio;
             cellnucratios(k) = cellnucheightratio;
             nucbottomslices(k) = nucbottomslice;
        else
            disp(['Illegal cell parameterization. Skipping ' fname])
        end
    else
        error([cellcodepath ' does not exist, cannot compute radius ratio.']);
    end
 end

% rad_ratio = distratio';
rad_ratio = distratio';

try
    %D. Sullivan 6/9/13 added support for flexible temp dir
    save([cellcodepath filesep 'radius_ratio.mat'],'rad_ratio','cellnucratios','nucbottomslices')

catch
end
