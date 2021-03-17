function answer = demo3D47()
% demo3D47
%
% Combine two generative model files into a single file.
%
% Input 
% -----
% * a list of valid CellOrganizer model files
%
% Output
% ------
% * a valid model

% Rita Chen
%
% Copyright (C) 2012-2018 Murphy Lab
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT MODIFY THIS BLOCK

if ~isdeployed()
  current_path = which(mfilename);
  [current_path, filename, extension] = fileparts( current_path );
  cd(current_path);
end

disp( 'demo3D47' );
disp( 'The estimated running time is about 1 minutes. Please wait.' );

options.output_filename = 'model.mat';
options.selection = [1,1,0;0,0,1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% documentation
% -------------
options.documentation.author = 'Murphy Lab';
options.documentation.email = 'murphy@cmu.edu';
options.documentation.website = 'murphy@cmu.edu';
options.documentation.description = 'This combined model is the result from demo3D47.';
options.documentation.date = date;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FEEL FREE TO MODIFY THE VARIABLES IN THIS BLOCK

% note, the framework will be synthesized using the first protein model
% found to match the given patterns in the SBML file.
% Changing the order/priority of this is not supported at this time.
directory = "../../../models/3D/";
file1 = strcat(directory,"mit.mat");
file2 = strcat(directory,"nuc.mat");
files = {file1, file2};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

answer = slml2slml(files,options);
end%demo3D47
