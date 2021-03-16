function protfit = protfit_percell(segdna, segcell, improt,immask,options)
%Trains protein patterns per cell and saves temporary files 

%Author: Devin Sullivan 6/7/13-6/13/13 - adapted from Tao Peng's 
%                                        train_protein_model2.m
%
%
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
% For additional information visit http://murphylab.web.cmu.edu/ or
% send email to murphy@cmu.edu


[objects, centers, blocksize, numBlocksX, numBlocksY, options] = findobjs_run_percell(improt,...
                                                                                    immask,...
                                                                                    options);

[mixes, objintens, offsets, options] = learngaussobjs_run3(objects, blocksize, options);

[objsize,intensities,objnum, param] = gmm_objempdistr_percell(mixes, objintens, options);
%             model.size = gmm_objsizefit(objstatdir);

size = gmm_objsizefit_percell(objsize);


tic;
%D. Sullivan 2/22/13 added param input to objpos_run to pass the
%resolution of the prot image
% model.position.beta = objpos_run(imgdir,savedir);
[beta, X, Y, imageID] =  objpos_run_percell(segdna, segcell, mixes, offsets, options);


protfit.objects = objects;
protfit.blocksize = blocksize;

gmm.mixes = mixes;
gmm.objintens = objintens;
gmm.offsets = offsets;

protfit.gmm = gmm;
protfit.size_model = size;
protfit.gauss_objsize = objsize;
protfit.gauss_intens = intensities;

protfit.pos.beta = beta;
protfit.pos.X = X;
protfit.pos.Y = Y;
toc

