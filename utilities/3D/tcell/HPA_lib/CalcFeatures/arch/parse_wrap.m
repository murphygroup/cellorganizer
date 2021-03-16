clear

addpath ./lib/system
addpath ./lib/fileio

fname = 'data/metadata/hpatable.mat';
if exist(fname,'file')
    load( fname)
else
    csvfilename = 'data/metadata/Murphy_ALL.csv';
    fid = fopen(csvfilename,'r');
    f = fread(fid);
    f = char(f)';
    fclose(fid);
    idx = find(f==13);
    f(idx) = [];

    mout = listmatlabformat( f);


    columns = parseline( mout{1}, ',');
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
    save(fname, 'hpatable','columns');
end

idx = find(hpatable.if_staining_intensity_id == 4 ...
    | hpatable.unspecific == 1 ...
    | hpatable.if_staining_intensity_id == 3 ...
    | hpatable.loc_cytoplasm == 1 ...
    | hpatable.stain_granular == 1);

for i=1:length(columns)
    hpatable.(columns{i})(idx) = [];
end

disp(['The number of samples that will be used is: ' ...
    num2str(length(hpatable.antibody))]);
% notes as of 10/16/2008: 6261 images, 5115 w/ filter, 1330 w/ stringent

abids = zeros(length(hpatable.antibody),1);
imglist = [];
imglist{length(hpatable.antibody)} = [];
rootdir = '/images/HPA/images/IFconfocal';
for i=1:length(hpatable.antibody)
    abids(i) = str2double( hpatable.antibody{i}(4:end));

    imglist{i} = [rootdir '/' ...
        num2str(hpatable.plate_id(i)) '/' ...
        num2str(hpatable.plate_id(i)) '_' ...
        hpatable.well{i} '_' ...
        '*_green.tif'];
end

[cellist I cellindex] = unique(hpatable.celline);



idx = strmatch( 'loc_', columns);
classes = columns(idx);
classmatrix = zeros(length(imglist),length(classes));
for i=1:length(classes)
    classmatrix(:,i) = hpatable.(classes{i});
end

% example for how to filter based on cell type
% classmatrix(cellindex==1,:);
% imglist(cellindex==1)';
% abids(cellindex==1);

% selecting only antibodies that stain one organelle
fname = 'data/metadata/hpalists.mat';
if exist(fname,'file')
    load(fname);
else
    specificity = sum(classmatrix,2);
    idx = find(specificity==1);

    mout = [];
    mout{length(idx)*3} = [];
    ind = zeros(length(idx)*3,1);
    st = 1;
    for i=1:length(idx)
        [uerr uout] = unix( ['ls -1 ' imglist{idx(i)}]);
        tout = listmatlabformat( uout);
        mout(st:st+length(tout)-1) = tout;
        ind(st:st+length(tout)-1) = idx(i);
        st = st+length(tout);
    end

    ind(st:end) = [];
    mout(st:end) = [];
    [a idx] = max(classmatrix(ind,:),[],2);

    imagelist = mout';
    antibodyids = abids(ind);
    cellabels = cellindex(ind);
    classlabels = classes(idx)';

    save(fname, 'imagelist','antibodyids','classlabels','cellabels');
end


% Watch out for dupe entries in imglist (only different ensid). ignoring
%    for now because img parsing/locking during processing should address this
