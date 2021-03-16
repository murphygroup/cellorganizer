function [img,imXYZ,resolution]  = model2microtubules( cellimg, nucimg, imgcent_coordinate, model, param )
% MICROTUBULES This function generates a microtubule from a given 
% model and parameters from a data set (datamat)
%
% Inputs:
% cell = cell shape model(negative, e.g. 0's inside cell)
% nuc = nuclear shape model (negative e.g. 0's inside nuc)
% imgcent_coordinate = coordinate
% datamat = matrix of microtubule parameter values from a given dataset
% savepath = where to save 
% microscope = name of microscope used to generate psf
% Outputs:
% image 

% Author: Devin Sullivan (devins@cmu.edu)
% Edited: Ivan E. Cao-Berg (icaoberg@scs.cmu.edu)
%
% March 8, 2012 Devin Sullivan Changed call to MT_synth
% March 18, 2012 Devin Sullivan floor number of microtubules to whole
% number
% Feb. 4, 2013 Devin Sullivan returned the resolution at which the image was
% generated (microns/pixel)
% Feb 9, 2013 Devin Sullivan returned imXYZ which was already being
% returned from MT_synth. This contains the XYZ coordinates of each MT as a
% cell array of MT's
% July 1, 2014 Ivan E. Cao-Berg Modified method so that it saves the microtubule
% temp results in the temp folder. This is very hacky and should be optimized before
% next release
%
% Copyright (C) 2012-2014 Murphy Lab
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

img = [];

cellimg = logical(~cellimg);
nucimg = logical(nucimg);
Dboth = cellimg+nucimg;

if nargin == 4
  microscope = 'none';
  XYZres = 0.2;
else
  try
    microscope = param.microscope;
  catch
    microscope = 'none';
  end

  try
    XYZres.cell = param.resolution.cell;
    XYZres.objects = model.resolution;
    XYZres.centrosome = param.resolution.centrosome;
  catch
    XYZres = 0.2;
  end
end 

%find the mean and covariance for each parameter
try
  meandata = model.parameters.mean;
catch
  warning('CellOrganizer: Incomplete model. Parameter mean not present or invalid.');
  return;
end

try   
  covmat = model.parameters.cov;
catch
  warning('CellOrganizer: Incomplete model. Parameter covariance not present or invalid.');
  return;
end

%take a random sample from the multivariate gaussian defined by the mean
%and covariance matrix from the data

%D. Sullivan 2/5/13 
%Need to look into how resolution effects these
sample = mvnrnd(meandata,covmat);
n = floor(sample(1));
mu_len = sample(2);
collin_min_number = sample(3);

%synthesize a microtubule distribution 
try
    %D. Sullivan 2/6/13 added resolution output to MT_synth
   [img,imXYZ,mtXYZ,randlengths,resolution] = MT_synth(n,mu_len,collin_min_number,Dboth,imgcent_coordinate,XYZres);
   
   %icaoberg 07/01/2014
   %modified code to save the temp results in the temp folder
   %this is very hacky and should be optimized before next release√
   temp_directory = [ pwd filesep 'temp' filesep 'microtubules' ];
   if ~exist( temp_directory )
        mkdir( temp_directory )
   end

   try
	index = length( dir( [temp_directory filesep '*.mat'] ) ) + 1;
   catch
	index = 1;
   end
	
   filename = [ temp_directory filesep 'microtubule' num2str(index) '.mat' ];

   %icaoberg
   save( filename, 'img', 'Dboth', 'imXYZ', 'mtXYZ', 'randlengths', 'resolution' );
catch err
   warning('CellOrganizer: Unable to synthesize a microtubule distribution.' );
   getReport( err )
   img = [];
   return;
end
