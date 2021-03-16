function tz_imshowdir_ds(drugroot)
%TZ_IMSHOWDIR_DS Show projection of 3D images for DS.
%   TZ_IMSHOWDIR_DS(DRUGROOT) shows the projection of 3D images under the
%   directory DRUGROOT, which should have the structure for saving drug
%   images.

%   ??-???-???? Initial write T. Zhao
%   31-OCT-2004 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

load diroot.mat

protlist=ml_dir([drugroot '/' DIDIR]);
protlist=protlist(3:end);
nprot=length(protlist);

for i=1:nprot
    protdir=[drugroot '/' DIDIR '/' protlist{i}];
    druglist=ml_dir(protdir);
    druglist=druglist(3:end);
    ndrug=length(druglist);
    
    for j=1:ndrug
        drugdir=[protdir '/' druglist{j}];
        celllist=ml_dir(drugdir);
        celllist=celllist(3:end);
        ncell=length(celllist);
        
        for k=1:ncell
            celldir=[drugdir '/' celllist{k}];
            maskpath=[drugroot '/' DMDIR '/' protlist{i} '/' druglist{j} '/' celllist{k} '.mat'];
            load(maskpath);
            tz_imshowdir(celldir,'tif','proj',double(mask));
        end
    end
end
