function stats = count_number_of_objects_wrapper( directory )

%COUNT_NUMBER_OF_OBJECTS Helper function that displays the number of
%objects per image as well as a some useful statistics from the
%intermediate results folder
% 
% Example
% >> directory = '~/demo3D12/temp/protein_objects_gaussian/original_objects';
% >> stats = count_number_of_objects( directory );

% Author: Ivan E. Cao-Berg (icaoberg@cmu.edu)
%
% Copyright (C) 2014 Murphy Lab
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

files = dir([directory filesep 'obj*.mat']);

names = {};
number_of_objects = [];

for index=1:1:length(files)
    file = files(index).name;
    
    try
        load( file );
    catch
        warning(['Unable to load file: ' file]);
    end
    
    try
    names{length(names)+1} = file;
    number_of_objects(length(number_of_objects)+1) = ...
        length(objects{1});
    disp([file, ': ', num2str(length(objects{1}))]);
    catch
        warning(['Loaded file: ' file ' ']);
    end
end

disp(['minimum: ' num2str(min(number_of_objects))]);
disp(['mean: ' num2str(mean(number_of_objects))]);
disp(['median: ' num2str(median(number_of_objects))]);
disp(['mode: ' num2str(mode(number_of_objects))]);
disp(['max: ' num2str(max(number_of_objects))]);

try
    hist( number_of_objects )
catch
    warning( 'Unable to make display' );
end

stats.names = names;
stats.number_of_objects = number_of_objects;
stats.minimum = min(number_of_objects);
stats.mean = mean(number_of_objects);
stats.median = median(number_of_objects);
stats.mode = mode(number_of_objects);
stats.max = max(number_of_objects);
