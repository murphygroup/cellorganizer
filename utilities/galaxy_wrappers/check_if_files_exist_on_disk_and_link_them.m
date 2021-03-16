function answer = check_if_files_exist_on_disk_and_link_them( string )
% CHECK_IF_FILES_EXIST_ON_DISK Simple helper function that takes a string
% from Galaxy and checks if these strings match local filenames

% Ivan E. Cao-Berg (icaoberg@cmu.edu), Xin Lu (xlu2@andrew.cmu.edu)
%
% Copyright (C) 2017 Murphy Lab
% Computational Biology Department
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

delimiter = ',';
C = strsplit(string,delimiter);

if isempty(C)
    answer = false;
    return
end

answer = true;
for i=1:length(C)
    filename = C{i};
    if exist(filename)
        disp(['File ' filename ' exists. Linking file to ' pwd filesep 'ome' num2str(i) '.tif']);
        system(['ln -s ' filename ' ' pwd filesep 'cell' num2str(sprintf('%05d',i)) '.ome.tif']);
    else
        warning(['File ' filename ' does not exist.']);
        answer = false;
    end
end
