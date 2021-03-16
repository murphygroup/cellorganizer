function answer = organize_synthesis_temp_files(index, options)
% ORGANIZE_SYNTHESIS_TEMP_FILES Helper function that movies and copies
% around temporary files to benefit debugging mode

% Copyright (C) 2016 Murphy Lab
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

answer = false;
output_directory = [options.temporary_results ...
    filesep 'cell' num2str(index) ];
if ~exist( output_directory )
    mkdir( output_directory );
end

files = dir([options.temporary_results filesep '*.mat']);
for i=1:1:length(files)
    filename = [options.temporary_results filesep files(i).name];
    if options.display
        disp( ['Moving temporary file ' filename] );
    end
    movefile( filename, output_directory );
end

files = dir([pwd filesep '*meshdata.mat']);
for i=1:1:length(files)
    filename = [pwd filesep files(i).name];
    if options.display
        disp( ['Moving temporary file ' filename] );
    end
    movefile( filename, output_directory );
end

answer = true;