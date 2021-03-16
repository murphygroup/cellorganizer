function [MTimage,MTimagewpsf] = microtubules(Dboth,imgcent_coordinate,microscope,datamat,XYZres)
%MICROTUBULES This function generates a microtubule from a given model and parameters from a data set (datamat)
%
% Inputs:
% Dboth = cell+nuc;
%   cell = cell shape model(negative, e.g. 0's inside cell)
%   nuc = nuclear shape model (negative e.g. 0's inside nuc)
% imgcent_coordinate = coordinate
% microscope = name of microscope used to generate psf
% datamat = matrix of microtubule parameter values from a given dataset
% XYZres = resolution of the pixels in microns
%
% Outputs:
% MTimage 
% MTimagewpsf
%
% note: this program assumes that the imgcent_coordinate is taken from the
% cell image and not the nuclear image and allows the cell and nuclear
% images to be potentially different sizes

% Author: Devin Sullivan (devins@cmu.edu)
%
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

%need to resize the cell and nuc images (downsample) for speed
%first make sure the image has the same x by y dimensions
if(size(Dboth,1)~=size(Dboth,2))
    diff = size(Dboth,1)-size(Dboth,2);
    disp('Cell image x and y dimension are not equal, adding cushion');
    if(diff>0)
        Dboth = [Dboth,zeros(size(Dboth,1),diff,size(Dboth,3))];
    elseif(diff<0)
        Dboth = [Dboth;zeros(diff,size(Dboth,2),size(Dboth,3))];
    end
end

%do the downsampling
if(size(Dboth,1)==256)
    disp('No downsampling required, image size is 256x256xZ');
else
   disp(['Downsampling image from ' size(cell) ' to 256x256x(Znew)']);
    Dboth = imresize(Dboth,(256/size(Dboth,1)));
    %not sure this is correct
    disp('Adjusting imgcent_coordinate');
    imgcent_coordinate = floor(imgcent_coordinate.*(256/size(Dboth,1)));    
end

if~exist('XYZres','var')
    disp('No resolution specified, using default value from HeLa dataset');
    XYZres = 0.2
end
%if no datamat is provided
if ~exist('datamat','var')
    disp('No data provided for microtubule parameters, using HeLa parameters');
    load('HeLaresult.mat','finTable');
    datamat = finTable(:,[1,2,4]);
end

%find the mean and covariance for each parameter
meandata = mean(datamat);
covmat = cov(datamat);

%take a random sample from the multivariate gaussian defined by the mean
%and covariance matrix from the data
sample = mvnrnd(meandata,covmat);
n = sample(1);
n = ceil(n);
mu_len = sample(2);
mu_len = ceil(mu_len);
collin_min_number = sample(3);
collin_min_number = ceil(collin_min_number*100)/100;


%sythesize a microtubule distribution 
disp('Generating microtubule image');
[MTimage,imXYZ,mtXYZ,randlengths] = MT_synth_dps(n,mu_len,collin_min_number,Dboth,imgcent_coordinate,XYZres);

disp('Fetching point spread function');
%get the psf and apply it
psffile = which_psf(microscope);
if(length(psffile)<1)
    disp('No point spread function specified or found, returning original');
else
  MTimagewpsf = psf_blur_hela_mean_dps(MTimage,psffile);    
end

%upsampling image to 1024x1024 for saving 
%all returned images should be 256x256xZ unless unsupported flag is thrown
    disp(['Resizing image for display: ' int2str(size(MTimage)),'to 1024x1024x(Znew)']);
    MTimage = imresize(MTimage,(1024/size(MTimage)));
    if(exist(MTimagewpsf,'var'))
        MTimagewpsf = imresize(MTimagewpsf,(1024/size(MTimage)));    
    end
    %shouldn't need to resize these, but just incase you want them to match
    %MTimg/eachother, e.g. if you originally passed a 512x512 cell/nuc. 
    Dboth = imresize(Dboth,4);
    
%%%%%
%optional, uncomment if you want it
%Show the image in a pretty surface plot
%%need to also define these two variables
%segcell = cell shape 
%segdna = nuc shape
%show_3D_image_dps(n,mu_len,colli_min_number,segdna,segcell,psffile,imagepath)
%%%%%

%save the generated images 
save(['outputs/images/new_syn_n_' num2str(n) '_mulen_' num2str(mu_len) '_colli_' num2str(collin_min_number) '.mat'], ...
      'MTimage','MTimagewpsf','microscope'); 
