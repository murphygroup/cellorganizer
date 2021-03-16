function [rad_ratio,nucdist,nucelldist] = extract_radius_ratio(cellcodes, equatorZ)
% From the cellcodes, extract the radius ratio of nuclear radius over cell
% radius under cylindrical systerm

% Author: Tao Peng
% Edited: Ivan E. Cao-Berg

% June 10 2013 D. Sullivan This method now precomputes radius ratio 
%                          statistics per-cell 
%                          extract_radius_ratiomodel compiles the model
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
k = 0;

% Calculate radius ratios and combining them
distratio = []; cellnucratios = []; nucbottomslices = [];
% celllist = ml_dir([cellcodepath filesep 'cellcodes_*']);

% for i = 1:length(celllist)
%     fname = [cellcodepath filesep celllist{i}];

k = k + 1;

baseplane(k) = equatorZ;
nucdist = [];
nucelldist = [];
stack = [];
for s = 1:length(cellcodes)
%            disp(['s:' num2str(s)]);
    if ~isempty(cellcodes{s})
        curr_nuc = cellcodes{s}.nucdist;
        curr_cell = cellcodes{s}.nucelldist;
        if length(curr_nuc) ~= 360
            h = -pi:(2*pi/length(curr_nuc)):pi;
            curr_nuc(end+1) = curr_nuc(1);
            curr_nuc = interp1(h,curr_nuc,H);
            curr_nuc(end) = [];
            curr_cell(end+1) = curr_cell(1);
            curr_cell = interp1(h,curr_cell,H);
            curr_cell(end) = [];
        end
        nucdist = [nucdist;curr_nuc];
        nucelldist = [nucelldist;curr_cell];
    end         
end
r = nucdist ./ nucelldist;
distratio = r(:);
rad_ratio = distratio';

end
