function load_wrap(setname)

addpath ./lib/segmentation/watershed/
addpath ./lib/segmentation/
addpath ./lib/process/
addpath ./lib/system/


ANTIBODY_IDS = {...
	'6669','447','992','5221','1523','8021','2141','7261',...
	'4497','3358','1677','248','4479','5861','4055','4334' ...
	};
CLASSLABELS = {...
        'nucleus','nucleolus','golgi','ER','mitochondria',...
          'endosome','lysosome','cytoskeleton',...
        'nucleus','nucleolus','golgi','ER','mitochondria',...
          'endosome','lysosome','cytoskeleton' ...
        };

features = [];
antilabels = [];
classlabels = [];
imagelabels = [];
clear imagelist;
counter = 0;
for zed = 1:length(ANTIBODY_IDS)

rootdir = ['../data/2_unmixedImages/' ANTIBODY_IDS{zed} '/'];
filetype = 'png';

uout = unixfind( rootdir, filetype);

readlist = listmatlabformat( uout);

mout = findreplacestring( uout, '/', '_');
mout = findreplacestring( mout, '.png', '.mat');
mout = findreplacestring( mout, '.._data_2_unmixedImages_', ...
    ['./data/features/' setname '/']);

readlist_features = listmatlabformat( mout);

f_tmp = [];
i_tmp = [];

for i=1:length(readlist_features)
    counter = counter+1;
    load( readlist_features{i});
    f_tmp = [f_tmp; feats];

    idx_ = find(readlist_features{i}=='_');
    imname = readlist_features{i}(idx_(end-2)+1:end-4);
    imagelist{counter} = imname;
    i_tmp = [i_tmp; counter*ones(size(feats,1),1)];
end
features = [features; f_tmp];
anti = str2num(ANTIBODY_IDS{zed});
antilabels = [antilabels; anti*ones(size(f_tmp,1),1)];
clear cltmp
cltmp(1:size(f_tmp,1)) = CLASSLABELS(zed);
classlabels = [classlabels cltmp];
imagelabels = [imagelabels; i_tmp];
end

classlabels = classlabels';

save(['data/results/features_' setname '.mat'] ,...
    'features','antilabels','classlabels','imagelabels','imagelist');

