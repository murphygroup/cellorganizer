function answer = demo3D57()
% demo3D57
%
% This demo illustrates using CellOrganizer to show protein enrichment plot
% for certain regions of the 3D T cell following the approach described in
%
% K. T. Roybal, T. E. Buck, X. Ruan, B. H. Cho, D. J. Clark, R. Ambler,
% H. M. Tunbridge, J. Zhang, P. Verkade, C. Wülfing, and R. F. Murphy (2016)
% Computational spatiotemporal analysis identifies WAVE2 and Cofilin as
% joint regulators of costimulation-mediated T cell actin dynamics.
% Science Signaling 9:rs3. doi: 10.1126/scisignal.aad4149.
%
%ShowTcellEnrichment
% Input
% -----
% * a set of t cell models with different time points
%
% Output
% ------
% * Plots of enrichment for different purposes.

% Author: Xiongtao Ruan
%
% Copyright (C) 2016-2019 Murphy Lab
% Computational Biology Department
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

if ~isdeployed()
  current_path = which(mfilename);
  [current_path, filename, extension] = fileparts( current_path );
  cd(current_path);
end
disp( 'demo3D57' );


filename = '../demo3D48/models';
if ~exist( filename )
    warning( [ 'File ' filename ' not found. Run demo3D48 first.' ] );
    answer = false;
    return
end

options = struct();
% Define the cylinder for the larger region (numerator)
% the radius of the ring
options.model.tcell.region_2_radius = [8, 8 / 3];
% the height of the cylinder
options.model.tcell.region_2_thickness = [4, 4];
% the top index of the cylinder starting from the top of the cell (1 means
% from the topmost of the cell)
options.model.tcell.region_2_start_ind = [1, 1];

% define the enrichment region as the ring
options.model.tcell.enrichment_region_type = 'ring';

% set the default parameter using top n% fluorescence as false
options.model.tcell.should_use_global_enrichment_region = false;

% set the paramter for user defined enrichment region as true
options.model.tcell.should_use_user_defined_enrichment_region = true;

% set the paramter for enrichment over certain region as true
options.model.tcell.enrichment_over_certain_region = true;

% set the parameter for the type of the region in the denominator
options.model.tcell.enrichment_bottom_region_type = 'top_fluorescence';

% set the errorbar type as sem
options.model.tcell.error_bar_type = 'sem';

% set the result filename for the all the enrichment
options.model.tcell.save_result_filename = 'enrichment_result_ring_over_top90.csv';

disp('Show enrichment plots for T cell model');
slml2info({'../demo3D48/models/LAT_reltime_-1.mat','../demo3D48/models/LAT_reltime_0.mat','../demo3D48/models/LAT_reltime_1.mat','../demo3D48/models/LAT_reltime_4.mat'},options)
end
