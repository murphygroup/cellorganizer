function frequency = getProteinFrequencyModel( savedir )
%GETPROTEINFREQUENCYMODEL Returns the protein frequency model component
%using the intermediate results calculated by CellOrganizer.

% Author: Ivan E. Cao-Berg (icaoberg@scs.cmu.edu)
%
% Edited: D. Sullivan 6/12/13 - refactored the code. now reading objects
%                               in from the object_stats folder. 
%                               (contains object intensities and sizes,
%                               + totnumber for each cell)
% Copyright (C) 2012 Murphy Lab
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

% April 24, 2012 R.F. Murphy Count gaussian objects found in each major
% object, not the major objects

%D. Sullivan 6/12/13 - refactored the code. now reading objects in from the
%object_stats folder. (contains object intensities and sizes, + totnumber)
% for j=1:intmax
    files = ml_dir( [savedir filesep 'gaussobjs_*.mat'] );
%     if isempty( files )
%       break
%     end
    numberOfObjects = zeros(1,length(files));
    for i=1:length(files)
        load( [savedir filesep files{i}],'objnum' );
        numberOfObjects(i) = objnum;
    end
% end

numberOfObjects=numberOfObjects(find(numberOfObjects~=0));

[parmhat,parmci] = lognfit( numberOfObjects );
frequency = [];
frequency.mu = parmhat(1);
frequency.sigma = parmhat(2);
