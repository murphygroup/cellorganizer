function ml_featcalc_demo
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


clc;
addpath(genpath('/home/ejr/ml/input/'));
addpath(genpath('/home/ejr/ml/featcalc/'));
defaultDir = '/imaging/mb_hela_images/mb_hela_images/giant/';
defaultBase = '/imaging/mb_hela_images/mb_hela_images/';
defaultIm = 'prot/r30jun97.giant.14--1---2.dat';
defaultDna = 'dna/r30jun97.gidap.14--1---2.dat';
defaultCrop = 'crop/30jun97.gicrop.14.tif';
defaultFeatSet = 'all84';

fprintf('  ----------------------------------------------------------\n');
fprintf('  |         DEMO for ml feature calculation code           |\n');
fprintf('  |         Press ENTER to use default arguments           |\n');
fprintf('  |           Default dir, image and feature set:          |\n');
fprintf('  |                                                        |\n');
fprintf('  |  %-54s|\n', defaultDir);
fprintf('  |      %-50s|\n', defaultIm);
fprintf('  |      %-50s|\n', defaultFeatSet);
fprintf('  |                                                        |\n');
fprintf('  ----------------------------------------------------------\n');


cont = 1;
currentFeatSet = 0;
currentImage = 0;
baseDir = defaultBase;
sets = [];
images = [];
while(cont)
    com = input('Command (h for help): ', 's');
    if (isempty(com)) continue; end
    switch com
     case {'h' 'H' 'help' 'Help'}
	fprintf('\t{q, Q, Quit, quit} - Quit\n');
	fprintf('\t{h, H, help, Help} - Help\n');
	fprintf(['\t{d, display, disp} - Display Features in ' ...
		 'Memory\n']);
	fprintf('\t{p, picts}         - Show Images in Memory\n');
	fprintf('\t{l, load}          - Load features\n');
	fprintf(['\t{f, feature}       - Show particular feature' ...
		 '\n']);
	fprintf('\t{Clear, clear, clc}- Clear Screen\n');
	fprintf('\t{c, comp}          - Compare two features\n');
	fprintf('\t{i, image}         - Load Image\n');
	fprintf('\t{s, show}          - Show Image\n');
	fprintf('\t{a, setBase}       - Set Base Directory\n');
	fprintf('\t{b, base}          - Show base directory\n');
	
	
     case {'q' 'Q' 'Quit' 'quit'}
	cont = 0;
     case {'d' 'display' 'disp'}
	if (length(sets) == 0)
	   fprintf('\tNo features have been calculated.\n');
	   continue;
	end
	fprintf('|'); 
	for i = 1:length(sets) 
	    fprintf('%10.0f|', i); 
	end
	fprintf('\n');
	fprintf('|');
	for i = 1:length(sets)
	    fprintf('%10.0f|', sets{i}.imageNum);
	end
	fprintf('\n');
	fprintf('|'); 
	for i = 1:length(sets) 
	    fprintf('%10s|', sets{i}.featSet); 
	end
	fprintf('\n');
	fprintf('Current set is %.0f\n', currentFeatSet);
     case {'a' 'setBase'}
	baseDir = input(sprintf('\tEnter Base directory: ', 's'));
     case {'b' 'base'}
	fprintf('%s\n', baseDir);
     case {'p' 'picts'}
	if (length(images) == 0)
	    fprintf('\tNo images have been loaded.\n');
	    continue;
	end
	fprintf('--------------------\n');
	for i = 1:length(images)
	    fprintf('%2.0f:\n', i);
	    fprintf('Image:\n%s\n', images{i}.im.name);
	    if (~isempty(images{i}.dna.name))
		fprintf('DNA:\n%s\n', images{i}.dna.name);
	    end
	    if (~isempty(images{i}.crop.name))
		fprintf('Crop:\n%s\n', images{i}.crop.name);
	    end
	    fprintf('--------------------\n');
	end
	fprintf('\n');
	fprintf('Current image is %.0f\n', currentImage);
     case {'l' 'load'}
	if (length(sets) > 5)
	    fprintf('\tSet limit reached.\n')
	    continue;
	end
	secondFeatSet = input(sprintf('\tEnter feature set name: '), ...
			      's');
	if (isempty(secondFeatSet) | secondFeatSet == ' ')
	    secondFeatSet = defaultFeatSet; 
	end
	[group] = calcFeats(images{currentImage}, secondFeatSet, currentImage);
	if (~isempty(group) & ~isa(group, 'char'))
	    sets{length(sets)+1} = group;
	    currentFeatSet = length(sets);
	end
     case {'f' 'feature'}
	featNum = input('Enter feature number: ', 's');
	featNum = str2num(featNum);
	displayFeat(featNum, images{currentImage});
     case { 'Clear' 'clear' 'clc'}
	clc;
     case { 'c' 'comp'}
	if (length(sets) < 2)
	    fprintf(['\tNot enough feature sets loaded for' ...
		     ' comparison!\n']);
	end
	one = input(sprintf('\tEnter first features(1-%2.0f): '), ... 
	                    length(sets), 's');
	two = input(sprintf('\tEnter second features(1-%2.0f): '), ...
		            length(sets), 's');

	if (one < 1 | one > length(sets) | two < 1 | two > ...
	    length(sets))
	   fprintf('\tInvalid feature selection!\n');
	   continue;
	end
	one = sets{one};
	two = sets{two};
	if (one.featSet ~= two.featSet)
	    fprintf(['\tInvalid feature selection! (Different' ...
		     ' features)\n']);
	end
	comp = input(sprintf('\tEnter feature to compare: '), 's');
     case {'i' 'image'}
	eval(['[group] = imageLoad([defaultDir defaultIm],' ...
	                          '[defaultDir defaultDna],' ...
	                          '[defaultDir defaultCrop], 1);'], ...
	     'fprintf(''Failed to load image.\n''); group = [];'); 
	if (~isempty(group))
	    images{length(images)+1} = group;
	    currentImage = length(images);
	end
     case {'s' 'show'}
	if (length(images)<1)
	    fprintf('\tNo images have been input.\n');
	    continue;
	end
	all = [];
	temp = images{currentImage};
	i = 1;
	all(:,:,1) = double(temp.im.image)/max(max(double(temp.im.image)));
	if (~isempty(temp.dna.image)) 
	    i = i+1;
	    all(:,:,i) = double(temp.dna.image)/ ...
		max(max(double(temp.dna.image)));
	end
	if (~isempty(temp.crop.image))
	    i = i+1;
	    all(:,:,i) = double(~temp.crop.image)/ ...
		max(max(double(~temp.crop.image)));
	end
	imshow(all);
     case {'demo'}
	
     otherwise
	if (com(1)== 'f')
	    featNum = str2num(com(2:end));
	    displayFeat(featNum, feat, name, slfname);
	elseif (com(1) == 'l')
	    if (length(sets) > 5)
		fprintf('\tSet limit reached.\n')
		continue;
	    end
	    secondFeatSet = com(2:end);
	    [group] = calcFeats(images{currentImage}, secondFeatSet, ...
				currentImage);
	    if (~isempty(group) & ~isa(group, 'char'))
		sets{length(sets)+1} = group;
		currentFeatSet = length(sets);
	    end
	elseif (com(1) == 's')
	    num = str2num(com(2:end));
	    if (isempty(num) | length(images)<num | num < 1)
		fprintf('%s is not a valid image.\n', com(2:end));
		continue;
	    end
	    if (length(images)<1)
		fprintf('\tNo images have been input.\n');
		continue;
	    end
	    all = [];
	    temp = images{num};
	    i = 1;
	    all(:,:,1) = double(temp.im.image)/max(max(double(temp.im.image)));
	    if (~isempty(temp.dna.image)) 
		i = i+1;
		all(:,:,i) = double(temp.dna.image)/ ...
		    max(max(double(temp.dna.image)));
	    end
	    if (~isempty(temp.crop.image))
		i = i+1;
		all(:,:,i) = double(~temp.crop.image)/ ...
		    max(max(double(~temp.crop.image)));
	    end
	    imshow(all);
        elseif (com(1) == 'x')
	    (eval(com(2:end), 'fprintf(''Invalid command: %s\n'', lasterr)'))
        end
	
    end
end

function [group] = imageLoad(imName, dnaName, cropName, def)
if def
    defaultIm = imName;
    defaultDna = dnaName;
    defaultCrop = cropName;
end


if (def | ~exist(imName))
    imName = input(sprintf('\tEnter image full path: '), 's');
    if (isempty(imName))
	imName = defaultIm;
    end
end

if (def | ~exist(dnaName))
    dnaCheck = input(sprintf('\tUse DNA image (y/n): ') , 's');
    if (isempty(dnaCheck))
	dnaName = defaultDna;
    else
	switch dnaCheck
	 case {'y' 'Y' 'Yes' 'yes' 't' 'T' 'True' 'true'}
	    dnaName = input(sprintf('\tEnter DNA full path: '), 's');
	    if (isempty(dnaName))
		dnaName = defaultDna;
	    end
	 otherwise
	    dnaName = [];
	end
    end
end

if (def | ~exist(cropName))
    cropCheck = input(sprintf('\tUse crop image (y/n): ') , 's');
    if (isempty(cropCheck))
	cropName = defaultCrop;
    else
	switch cropCheck
	 case {'y' 'Y' 'Yes' 'yes' 't' 'T' 'True' 'true'}
	    cropName = input('\tEnter crop full path: ', 's');
	    if (isempty(cropName))
		cropName = [defaultDir defaultCrop];
	    end
	 otherwise
	    cropName = [];
	end
    end
end


fprintf('Reading image:\n''%s''\n', imName);
im = ml_readimage(imName);
if (isempty(dnaName))
    dna = [];
else
    eval('dna = ml_readimage(dnaName);fprintf(''\twith DNA\n'');', ...
	 'fprintf(''\tIgnoring DNA\n''); dna = [];');
    
end

if (isempty(cropName))
    crop = [];
else
    eval('crop = ml_readimage(cropName);fprintf(''\twith cropping\n'');', ...
	 'fprintf(''\tIgnoring crop\n''); crop = [];');
end

group =  struct('im', struct('image', im, 'name', imName), ...
	        'dna', struct('image', dna, 'name', dnaName), ...
		'crop', struct('image', crop, 'name', cropName));



function displayFeat(featNum, feat, name, slfname)
if (isempty(featNum) | featNum > length(slfname) | featNum < 1)
    fprintf('\tInvalid feature\n');
else
    fprintf('\t%s(%s) = %d\n', name{featNum}, slfname{featNum}, ...
	    double(feat(featNum)));
end

function [group] = calcFeats(inIm, featSet, imNum)
im = inIm.im.image;
dna = inIm.dna.image;
crop = inIm.crop.image;
name = {};
feat = [];
slfname = {};
fail = 0;
fprintf('Running %s calculations on image.\n', featSet);
runString = ['[name, feat, slfname] = ml_featset(im, crop, dna,' ...
	     ' featSet);'];
errorString = ['fprintf(''Failed ml_featset: %s\n'', lasterr);' ...
	       'fail = 1;'];
t0 = cputime;
eval(runString, errorString);
t = cputime - t0;
if (~fail)
    fprintf('\n%i features calculated in approximately %.0f cpu seconds\n', ...
	    length(slfname), t);
    group = struct('name', {name}, 'feat', {feat}, 'slfname', ...
	           {slfname}, 'featSet', featSet, 'imageNum', imNum);
else
    fprintf('\nFeatures failed to calculate after %.0f cpu seconds\n', ...
	    t);
    group = 'FAIL';
end

