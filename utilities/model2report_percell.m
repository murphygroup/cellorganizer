function model2report_percell(param)
%Function that creates useful reports on percell statistics 
%Sometimes we may wish to use percell for comparision as in a movie 
%
%Inputs: 
%param = structure of params including temporary folder where the various
%        models are saved.
%
%Outputs: 
%
%A set of graphs and pvalues 

%Author: Devin Sullivan 6/17/13
%
% Copyright (C) 2007-2013  Murphy Lab
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

%%
%Objects 
%Figure out how many cell images there are 
files = ml_dir([param.objstatsdir filesep '*.mat']);


%initialize stats
meanint = zeros(1,length(files));
stdint = zeros(1,length(files));
totnum = zeros(1,length(files));
meansize = zeros(length(files),3);
stdsize = zeros(length(files),3);
meanvolume = zeros(1,length(files));
stdvolume = zeros(1,length(files));

for i = 1:length(meanint)
    load([param.objstatsdir filesep files{i}]);
    meanint(i) = mean(intensities);
    stdint(i) = std(intensities);
    meansize(i,:) = mean(objsize);
    stdsize(i,:) = std(objsize);
    totnum(i) = objnum;
    %volume
    vol = (4/3*pi).*(objsize(:,1)/.2.*objsize(:,2)./2.*objsize(:,3)./2);
    meanvolume(i) = mean(vol);
    stdvolume(i) = std(vol);
    
end

%adjust volume if you know the resolution to have a measurment in um
%instead of pixels
if isfield(param.model,'resolution');
    %convert the pixel volumes to micrometers 
    if size(param.model.resolution,2)==2 || size(param.model.resolution,2)==3
        meanvolume = meanvolume.*prod(param.model.resolution);
        ylabelstr = 'mean volume (um)';
    else
        disp(['Resolution specified but not recognized. Please input a '...
             '1x2 or 1x3 vector']);
    end
else
    ylabelstr = 'mean volume (pixels^3) - careful of non-cubic pixels';
end

%volume plot
try
    titlestr = 'Object volume';
    xlabelstr = 'Cell num';
    niceplot([1:length(files)],meanvolume,xlabelstr,ylabelstr,titlestr,'plot','.')
catch
    warning('Unable to open display' );
end

%number plot
try
    xlabelstr = 'Cell num';
    ylabelstr = 'number of objects/cell';
    titlestr = 'Number of objects';
    niceplot([1:length(files)],totnum,xlabelstr,ylabelstr,titlestr,'plot','.')
catch
    warning('Unable to open display' );
end
% figure,plot(1:length(files),totnum,'.','MarkerSize',20);
% title('Number of objects','FontSize',20)
% xlabel('Cell num','FontSize',20);
% ylabel('number of objects/cell','FontSize',20);
%%

%%
%Compartments 
files = ml_dir([param.compartmentdir filesep '*.mat']);

for i = 1:length(files)
    load([param.compartmentdir filesep files{i}]);
    propfluortot(i,:) = propfluor;
    
end

%proportional fluorescence plot 
xlabelstr = 'Cell num';
ylabelstr = 'Proportion of fluorescence';
titlstr = 'Compartmental fluorescent proportions';

try
    niceplot([1:length(files)],propfluortot,xlabelstr,ylabelstr,titlestr,'semilogy','.')
    legend(compartmentlist,'FontSize',20,'Fontname','Ariel')
catch
    warning('Unable to open display' );
end
