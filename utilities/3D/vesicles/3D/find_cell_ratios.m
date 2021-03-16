function allratios = find_cell_ratios( impath, savepath, param )
%FIND_CELL_RATIOS

% March 28, 2012 R.F. Murphy from find_cell_codes by T. Zhao
%
% Copyright (C) 2012 Murphy Lab
% Lane Center for Computational Biology
% School of Computer Science
% Carnegie Mellon University
%
% March 12, 2012 Added another preprocessing routine. This routine mimics
% the preprocessing from A. Shariff, G. K. Rohde and R. F. Murphy (2010) A
% Generative Model of Microtubule Distributions, and Indirect Estimation of
% its Parameters from Fluorescence Microscopy Images. Cytometry Part A 77A:457-466.
% July 26, 2012 I. Cao-Berg If the number of slices in the nuclear image that have
%                           fluorescence is less than 4, then ignore such image
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

if nargin == 2
   param = [];
end

try
   downsample = param.downsample;
catch
   downsample = [5 5 1];
end

celldirlist = ml_dir([impath '/cell*']);
tmpdir = [ pwd filesep 'temp' filesep 'preprocessing' ];
if ~exist( tmpdir, 'file' )
   mkdir( tmpdir );
end

allratios = [];
for i = 1:length(celldirlist)    
    savefile = [savepath  filesep 'cellratios_' int2str(i) '.mat'];
    if ~exist( savefile , 'file' )
        disp([ 'Image ' num2str(i)]);
        tic        
        try
            cellfile = [ tmpdir filesep 'cell' celldirlist{i}(5:end) '.mat' ];
            if ~exist( cellfile, 'file' )
               [segdna,segcell] = seg_cell_and_nucleus([impath filesep ...
                   celldirlist{i}], downsample );
               save( cellfile, 'segdna', 'segcell' );
            else
               load( cellfile );
            end

            %icaoberg 26/7/2012
            nnzpixel = squeeze(sum(sum(segdna)));

            if length(find(nnzpixel)>0) < 4
               continue
            end
            
            nucstacknumber = 11;
            [segdna, segcell] = ml_rescaleImage2Nucleus( segdna, segcell, nucstacknumber );
            ratios=ml_calccellnucratios3d(segcell,segdna,1);
            save([savepath filesep 'cellratios_' int2str(i)],'ratios');

            allratios = [allratios;ratios];
            toc;
        catch
            disp( ['Ignoring image ' num2str(i)] );
            toc;
        end
    else
       disp('Found and loading');
       load( [savepath filesep 'cellratios_' int2str(i) '.mat'] );
       allratios = [allratios;ratios];
    end
end
end%find_cell_ratios
