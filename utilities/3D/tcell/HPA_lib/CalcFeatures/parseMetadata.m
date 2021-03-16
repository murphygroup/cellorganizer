function parseMetadata(rootmeta,metafile,rootdir)

% rootmeta is the root directory storing metadata.
% metafile is the metadata file.
% rootdir is the directory storing all images.

% Copyright (C) 2010  Murphy Lab
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

% 10 Jan 10 - Jieyue Li


if ~exist('rootmeta','var')
    rootmeta = '/home/jieyuel/images/HPA/metadata/';
end

if ~exist('metafile','var')
    metafile = 'IF_atlas_5.csv';
end

if ~exist('rootdir','var')
    rootdir = '/images/HPA/images/IFconfocal';
end

fname = [rootmeta,'hpatable.mat']; %%
if exist(fname,'file')
    load( fname)
else
    csvfilename = [rootmeta,metafile]; %%
    fid = fopen(csvfilename,'r');
    f = fread(fid);
    f = char(f)';
    fclose(fid);
    idx = find(f==13);
    f(idx) = [];

    mout = listmatlabformat( f);


    columns = parseline( mout{1}, ',');
    columns = regexprep(columns,' ','_'); %%
    columns = regexprep(columns,'(',''); %%
    columns = regexprep(columns,')',''); %%
    %keyboard
    template = [];
    for i=1:length(columns)
        template = setfield(template,columns{i},[]);
    end

    hpatable = [];
    hparow = parseline( mout{2}, ',');
    columntypes = zeros(1,length(columns));
    for i=1:length(columns)
        columntypes(i) = str2double(hparow{i});
        if isnan(columntypes(i))
            array = [];
            array{length(mout)-1} = [];
            array = array';
        else
            array = zeros(length(mout)-1,1);
        end
        hpatable.(columns{i}) = array;
    end

    for i=2:length(mout)
        hparow = parseline( mout{i}, ',');
        for j=1:length(hparow)
            entry = hparow{j};
            if isnan(columntypes(j))
                hpatable.(columns{j}){i-1} = entry;
            else
                hpatable.(columns{j})(i-1) = str2double(entry);
            end
        end
    end
    
    % some repeat entries in the csv file, fixed 3/16/09
    str = [];
    for i=1:length(hpatable.Position)
        pid = num2str(hpatable.Plate_id(i));
        wid = hpatable.Position{i};
        str{i} = [rootdir pid '/' pid '_' wid '_*_green.tif'];
    end
    [u ind j] = unique(str);

    fn = fieldnames(hpatable);
    for i=1:length(fn)
        hpatable.(fn{i}) = hpatable.(fn{i})(ind);
    end
    
    save(fname, 'hpatable','columns');
end

%idx = find(hpatable.Staining_intensity == 4 ...
%    | hpatable.unspecific == 1);
idx = strcmp(hpatable.Staining_intensity,'4:Negative');

% ...
%    | hpatable.Staining_intensity == 3 ...
%    | hpatable.loc_cytoplasm == 1 ...
%    | hpatable.stain_granular == 1);

for i=1:length(columns)
    hpatable.(columns{i})(idx) = [];
end

disp(['The number of samples that will be used is: ' ...
    num2str(length(hpatable.Antibody))]); %%
% notes as of 10/16/2008: 6261 images, 5115 w/ filter, 1330 w/ stringent

abids = zeros(length(hpatable.Antibody),1);
imglist1 = [];
imglist1{length(hpatable.Antibody)} = [];
imglist2 = [];
imglist2{length(hpatable.Antibody)} = [];
imglist3 = [];
imglist3{length(hpatable.Antibody)} = [];

rootdir1 = rootdir; %%
rootdir2 = rootdir; %%
rootdir3 = rootdir; %%
for i=1:length(hpatable.Antibody)
    abids(i) = str2double( hpatable.Antibody{i}(4:end));

    imglist1{i} = [rootdir1 '/' ...
        num2str(hpatable.Plate_id(i)) '/' ...
        num2str(hpatable.Plate_id(i)) '_' ...
        hpatable.Position{i} '_' ...
        '*_green.tif'];
    imglist2{i} = [rootdir2 '/' ...
        num2str(hpatable.Plate_id(i)) '/' ...
        num2str(hpatable.Plate_id(i)) '_' ...
        hpatable.Position{i} '_' ...
        '*_green.tif'];
    imglist3{i} = [rootdir3 '/' ...
        num2str(hpatable.Plate_id(i)) '/' ...
        num2str(hpatable.Plate_id(i)) '_' ...
        hpatable.Position{i} '_' ...
        '*_green.tif'];
end

[cellist I cellindex] = unique(hpatable.Celline);



idx = strmatch( 'loc_', columns);
classes = columns(idx);
classmatrix = zeros(length(imglist1),length(classes));
for i=1:length(classes)
    classmatrix(:,i) = hpatable.(classes{i});
end


classes2 = [];
for i=1:size(classmatrix,1)
    cls = classes(find(classmatrix(i,:)));
    classes2{i} = '';
    for j=1:length(cls)
        classes2{i} = [classes2{i} cls{j} ','];
    end
    if length(classes2{i})
        classes2{i}(end) = [];
    end
end


%fname = '/home/jieyuel/images/HPA/metadata/hpalistsall.mat';
fname = regexprep(fname, 'hpatable.mat', 'hpalistsall.mat');
if exist(fname,'file')
    load(fname);
else
    specificity = sum(classmatrix,2);
    idx = 1:1:size(classmatrix,1);

    mout = [];
    mout{length(idx)*3} = [];
    ind = zeros(length(idx)*3,1);
    st = 1;
    for i=1:length(idx)
        [uerr uout] = unix( ['ls -1 ' imglist1{idx(i)}]);
        if uerr~=0
           [uerr uout] = unix( ['ls -1 ' imglist2{idx(i)}]);
           if uerr~=0
              [uerr uout] = unix( ['ls -1 ' imglist3{idx(i)}]);
           end
        end
        if uerr~=0
           i
           warning(['no such image ',imglist1{idx(i)},' or ',imglist2{idx(i)},' or ',imglist3{idx(i)}]);
           continue;
        end
        tout = listmatlabformat( uout);
        mout(st:st+length(tout)-1) = tout;
        ind(st:st+length(tout)-1) = idx(i);
        st = st+length(tout);
    end

    ind(st:end) = [];
    mout(st:end) = [];
    [a idx] = max(classmatrix(ind,:),[],2);

    imagelist = mout';
    %keyboard
    Antibodyids = abids(ind);
    cellabels = cellindex(ind);
    classlabels = classes2(ind)';
    specificity = specificity(ind);
    staining = hpatable.Staining_intensity(ind);

    antibodyids = Antibodyids;
    clear Antibodyids;

    save(fname, 'imagelist','antibodyids','classlabels','cellabels','specificity','staining');
end


% Watch out for dupe entries in imglist (only different ensid). ignoring
%    for now because img parsing/locking during processing should address this
