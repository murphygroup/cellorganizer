
addpath ./lib/system/
addpath ./lib/features/
addpath ./lib/segmentation/output/

rootdir = '/images/HPA/confocal/';
filetype = 'tif';
greparg = 'green';

ind = find(rootdir=='/');
rootdir_ = rootdir;
rootdir_(ind) = '_';

uout = unixfind( rootdir, filetype, greparg);
readlist = listmatlabformat( uout);

uout_nuc = findreplacestring( uout, 'green','blue');
uout_tub = findreplacestring( uout, 'green','red');
uout_er = findreplacestring( uout, 'green','yellow');
readlist_nuc = listmatlabformat( uout_nuc);
readlist_tub = listmatlabformat( uout_tub);
readlist_er = listmatlabformat( uout_er);


mout = findreplacestring( uout, '/', '_');
mout = findreplacestring( mout, '.tif', '.png');
mout = findreplacestring( mout, rootdir_, ['./data/masks/']);
readlist_mask = listmatlabformat( mout);

mout = findreplacestring( uout, '/', '_');
mout = findreplacestring( mout, '.tif', '.jpg');
mout = findreplacestring( mout, rootdir_, ['./data/maskvis/']);
writelist = listmatlabformat( mout);

for i=1:length(readlist)
    if exist(writelist{i},'file')
        continue;
    end

    tmpfile = writelist{i};
    tmpfile(find(tmpfile=='/')) = [];
    tmpfile(find(tmpfile=='.')) = [];
    tmpfile = ['./tmp/' tmpfile '.txt'];

    if exist(tmpfile,'file')
        continue;
    end
    fid = fopen(tmpfile,'w');

    nuc = imread(readlist_nuc{i});
    cell = imread(readlist_er{i});
    masks = imread(readlist_mask{i});
    output = segmentedfield(nuc,cell,masks);

    imwrite(output, writelist{i});

    fclose(fid);
    delete(tmpfile);
end

