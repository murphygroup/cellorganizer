function [combz,comby,combx]=ml_3dorthoview( imgs, viewtype)
%ML_3DORTHOVIEW Get three orthogonal slices of a 3D image with 3 channels
%   COMBZ = ML_3DORTHOVIEW(IMGS,VIEWTYPE) returns the slice with at the
%   center of fluorescence along Z direction if VIEWTYPE is 'sect'. It
%   returns maximum projection if VIEWTYPE is 'proj'. IMGS is a three
%   dimensional matrix. IMGS(:,:,1) corresponds to protein channel,
%   IMGS(:,:,2) corresponds to dna channel and IMGS(:,:,3) corresponds to
%   cell channel.
%   
%   [COMBZ,COMBY,COMBX] = ML_3DORTHOVIEW(...) also returns images along Y
%   and X directions.
%   
%   See also

% Copyright (C) 2006  Murphy Lab
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

% Background subtract, Crop and mask the images
%DNA
dnaCOF = ml_findCOF( ml_sparse( imgs{2}));
%dnathresh = 255*mb_nihthreshold( dnaclean);
%dnabin = ml_binarize( dnaclean, uint8(dnathresh));
%clear dnaclean;

%Cell
cellCOF = ml_findCOF( ml_sparse( imgs{3}));
%cellthresh = 15; %255*ml_threshold( cellclean);
%cellbin = ml_binarize( cellclean, uint8(cellthresh));
%clear cellclean;

%Prot
COF = ml_findCOF( ml_sparse( imgs{1}));
%protthresh = 255*ml_threshold( protclean);
%protbin = ml_binarize( protclean, uint8(protthresh));

protclean = imgs{1};
dnaclean = imgs{2};
cellclean = imgs{3};

cof = COF;
MaxVal = 1;

switch viewtype
    case 'sect'
        %     f1
        z=round(cof(3));
        slicez = stretch(double(protclean(:,:,z)),MaxVal);
        cellslicez = stretch(double(cellclean(:,:,z)),MaxVal);
        dnaslicez = stretch(double(dnaclean(:,:,z)),MaxVal);
        combz = cat(3,dnaslicez,slicez,cellslicez);
        imagesc(combz)
        %     gr
        truesize
        %     f2
        x = round(cof(2));
        slicex = stretch(double(permute(protclean(:,x,:),[1 3 2])),MaxVal);
        cellslicex = stretch(double(permute(cellclean(:,x,:),[1 3 2])),MaxVal);
        dnaslicex = stretch(double(permute(dnaclean(:,x,:),[1 3 2])),MaxVal);
        combx = cat(3,dnaslicex,slicex,cellslicex);
        imagesc(combx)
        %     gr
        truesize
        %     f3
        y = round(cof(1));
        slicey = stretch(double(permute(protclean(y,:,:),[3 2 1])),MaxVal);
        cellslicey = stretch(double(permute(cellclean(y,:,:),[3 2 1])),MaxVal);
        dnaslicey = stretch(double(permute(dnaclean(y,:,:),[3 2 1])),MaxVal);
        comby = cat(3,dnaslicey,slicey,cellslicey);
        imagesc(comby);
        %     gr
        truesize
        %r =input('h');
    case 'proj'
        %     f1
        projz = stretch(double(sum(protclean,3)),MaxVal);
        cellprojz = stretch(double(sum(cellclean,3)),MaxVal);
        dnaprojz = stretch(double(sum(dnaclean,3)),MaxVal);
        combz = cat(3,dnaprojz,projz,cellprojz);
        imagesc(combz);
        %     gr
        truesize
        %     f2
        projx = stretch(double(permute(sum(protclean,2),[1 3 2])),MaxVal);
        cellprojx = stretch(double(permute(sum(cellclean,2),[1 3 2])),MaxVal);
        dnaprojx = stretch(double(permute(sum(dnaclean,2),[1 3 2])),MaxVal);
        combx = cat(3,dnaprojx,projx,cellprojx);
        imagesc(combx);
        %     gr
        truesize
        %     f3
        projy = stretch(double(permute(sum(protclean,1),[3 2 1])),MaxVal);
        cellprojy = stretch(double(permute(sum(cellclean,1),[3 2 1])),MaxVal);
        dnaprojy = stretch(double(permute(sum(dnaclean,1),[3 2 1])),MaxVal);
        comby = cat(3,dnaprojy,projy,cellprojy);
        imagesc(comby);
        %     gr
        truesize
end
