function fieldFeatCalc

addpath lib/general/helper
addpath ./lib/system/
addpath ./lib/features/
addpath ./lib/segmentation/output/
addpath /home/jnewberg/CVS/SLIC
setcvspath('/home/jnewberg/CVS/SLIC/')

fsetnames = {...
    'texture16Field', ...
    'texture32Field', ...
    'texture128Field', ...
    'texture256Field', ...
    'texture16x2Field', ...
    'texture32x2Field', ...
    'texture128x2Field', ...
    'texture256x2Field', ...
    'texture16x4Field', ...
    'texture32x4Field', ...
    'texture64x4Field', ...
    'texture128x4Field', ...
    'texture256x4Field' ...
};

%     'skelField'... takes toooooo long
% fsetnames = {'objField'};


rootdir = '/images/HPA/images/IFconfocal/';
filetype = 'tif';
greparg = ' | grep green';

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

DSF = 2;
for i=1:length(readlist)
    prot = [];
    nuc = [];
    tub = [];
    er = [];
    prot_dsX = [];
    nuc_dsX = [];
    tub_dsX = [];
    er_dsX = [];
    thr_prot_dsX = [];
    fg_prot_dsX = [];
    thr_nuc_dsX = [];
    fg_nuc_dsX = [];
    thr_tub_dsX = [];
    fg_tub_dsX = [];
    thr_er_dsX = [];
    fg_er_dsX = [];
    prot_ds2 = [];
    prot_ds4 = [];
    for zed = 1:length(fsetnames)
        mout = findreplacestring( readlist{i}, '/', '_');
        mout = findreplacestring( mout, '.tif', ['_' fsetnames{zed} '.mat']);
        writepath = findreplacestring( mout, rootdir_, ['./data/features/field/']);
        
        if exist(writepath,'file')
            continue;
        end

        tmpfile = writepath;
        tmpfile(find(tmpfile=='/')) = [];
        tmpfile(find(tmpfile=='.')) = [];
        tmpfile = ['./tmp/' tmpfile '.txt'];

        if exist(tmpfile,'file')
            continue;
        end
        fid = fopen(tmpfile,'w');


        switch fsetnames{zed}
          case 'tasField'
            prot = loadChannel( readlist{i}, prot);
            [names, feats, slfnames] = ml_tas( prot,0);

          case 'tasx2Field'
            prot = loadChannel( readlist{i}, prot);
            prot_dsX = downsampleChannel( prot, 2, prot_dsX);
            [names, feats, slfnames] = ml_tas( prot_dsX,0);
            
          case {'mutualInfo16Field','mutualInfo32Field','mutualInfo64Field'}
            prot = loadChannel( readlist{i}, prot);
            nuc = loadChannel( readlist_nuc{i}, nuc);
            tub = loadChannel( readlist_tub{i}, tub);
            er = loadChannel( readlist_er{i}, er);
            
            NBINS = str2double( fsetnames{zed}(end-6:end-5));
            NBINS(2) = NBINS(1);

            sprot = size(prot);
            prot_re = reshape(prot,[sprot(1)*sprot(2) 1]);
            nuc_re = reshape(nuc,[sprot(1)*sprot(2) 1]);
            tub_re = reshape(tub,[sprot(1)*sprot(2) 1]);
            er_re = reshape(er,[sprot(1)*sprot(2) 1]);

            feats1 = mutualinformation( [prot_re nuc_re], NBINS);
            feats2 = mutualinformation( [prot_re tub_re], NBINS);
            feats3 = mutualinformation( [prot_re er_re], NBINS);
            feats4 = mutualinformation( [tub_re er_re], NBINS);

            feats1(1,2) = corr2( prot_re, nuc_re);
            feats2(1,2) = corr2( prot_re, tub_re);
            feats3(1,2) = corr2( prot_re, er_re);
            feats4(1,2) = corr2( tub_re, er_re);

            feats = [feats1 feats2 feats3 feats4];
            str = [ '_' num2str(NBINS(1)) '_graylevels'];
            names = {['nuc:mutual_information' str],['nuc:correlation' str], ...
                ['tub:mutual_information' str],['tub:correlation' str], ...
                ['er:mutual_information' str],['er:correlation' str], ...
                ['tub_er:mutual_information' str],['tub_er:correlation' str]};
            slfnames = repmat({''},[1 length(names)]);
            
          case {'texture16Field','texture32Field','texture64Field', ...
                'texture128Field','texture256Field'}
            prot = loadChannel( readlist{i}, prot);
            GLEVELS = str2double( fsetnames{zed}(end-7:end-5));
            if isnan(GLEVELS)
                GLEVELS = str2double( fsetnames{zed}(end-6:end-5));
            end
            [feats,names,slfnames] = features_texture( prot, GLEVELS);
            
          case {'texture16x2Field','texture32x2Field','texture64x2Field', ...
                'texture128x2Field','texture256x2Field'}
            prot = loadChannel( readlist{i}, prot);
            prot_ds2 = downsampleChannel( prot, 2, prot_ds2);
            GLEVELS = str2double( fsetnames{zed}(end-9:end-7));
            if isnan(GLEVELS)
                GLEVELS = str2double( fsetnames{zed}(end-8:end-7));
            end
            [feats,names,slfnames] = features_texture( prot_ds2, GLEVELS);
            
            for m=1:length(names)
                names{m} = [names{m} '_2x_downsampled'];
            end

          case {'texture16x4Field','texture32x4Field','texture64x4Field', ...
                'texture128x4Field','texture256x4Field'}
            prot = loadChannel( readlist{i}, prot);
            prot_ds4 = downsampleChannel( prot, 4, prot_ds4);
            GLEVELS = str2double( fsetnames{zed}(end-9:end-7));
            if isnan(GLEVELS)
                GLEVELS = str2double( fsetnames{zed}(end-8:end-7));
            end
            prot_ds2 = downsampleChannel( prot, 4, prot_ds4);
            [feats,names,slfnames] = features_texture( prot_ds4, GLEVELS);
            
            for m=1:length(names)
                names{m} = [names{m} '_4x_downsampled'];
            end
            
          case 'overlapField'
            prot = loadChannel( readlist{i}, prot);
            nuc = loadChannel( readlist_nuc{i}, nuc);
            tub = loadChannel( readlist_tub{i}, tub);
            er = loadChannel( readlist_er{i}, er);
            
            [fg_prot thr_prot] = thresholdChannel( prot, fg_prot, thr_prot);

            [feats1,names1] = features_overlap( prot, fg_prot, fg_nuc, 'nuc');
            [feats2,names2] = features_overlap( prot, fg_prot, fg_tub, 'tub');
            [feats3,names3] = features_overlap( prot, fg_prot, fg_er, 'er');

            feats = [feats1 feats2 feats3];
            names = [names1 names2 names3];
            slfnames = repmat({''},[1 length(names)]);
            
            case 'overlapFieldx2'
            prot = loadChannel( readlist{i}, prot);
            nuc = loadChannel( readlist_nuc{i}, nuc);
            tub = loadChannel( readlist_tub{i}, tub);
            er = loadChannel( readlist_er{i}, er);

            prot_dsX = downsampleChannel( prot, DSF, prot_dsX);
            nuc_dsX = downsampleChannel( nuc, DSF, nuc_dsX);
            tub_dsX = downsampleChannel( tub, DSF, tub_dsX);
            er_dsX = downsampleChannel( er, DSF, er_dsX);
            
            [fg_prot_dsX thr_prot_dsX] = thresholdChannel( prot_dsX, fg_prot_dsX, thr_prot_dsX);

            [feats1,names1] = features_overlap( prot_dsX, fg_prot_dsX, fg_nuc_dsX, 'nuc');
            [feats2,names2] = features_overlap( prot_dsX, fg_prot_dsX, fg_tub_dsX, 'tub');
            [feats3,names3] = features_overlap( prot_dsX, fg_prot_dsX, fg_er_dsX, 'er');

            feats = [feats1 feats2 feats3];
            names = [names1 names2 names3];
            slfnames = repmat({''},[1 length(names)]);


          case 'nonObjFluorField'
            prot = loadChannel( readlist{i}, prot);
            nuc = loadChannel( readlist_nuc{i}, nuc);
            tub = loadChannel( readlist_tub{i}, tub);
            er = loadChannel( readlist_er{i}, er);
            
            prot_dsX = downsampleChannel( prot, DSF, prot_dsX);
            nuc_dsX = downsampleChannel( nuc, DSF, nuc_dsX);
            tub_dsX = downsampleChannel( tub, DSF, tub_dsX);
            er_dsX = downsampleChannel( er, DSF, er_dsX);
            
            [fg_prot_dsX thr_prot_dsX] = thresholdChannel( prot_dsX, fg_prot_dsX, thr_prot_dsX);
            [fg_nuc_dsX thr_nuc_dsX] = thresholdChannel( nuc_dsX, fg_nuc_dsX, thr_nuc_dsX);
            [fg_tub_dsX thr_tub_dsX] = thresholdChannel( tub_dsX, fg_tub_dsX, thr_tub_dsX);
            [fg_er_dsX thr_er_dsX] = thresholdChannel( er_dsX, fg_er_dsX, thr_er_dsX);

            prot_int = sum(fg_prot_dsX(:));
            bg_prot_dsX = prot_dsX;
            bg_prot_dsX(prot_dsX>thr_prot_dsX) = 0;
            bg_int = sum(bg_prot_dsX(:));
            bg_prot_fg_nuc_int = bg_prot_dsX;
            bg_prot_fg_nuc_int(fg_nuc_dsX==0) = 0;
            bg_prot_fg_tub_int = bg_prot_dsX;
            bg_prot_fg_tub_int(fg_tub_dsX==0) = 0;
            bg_prot_fg_er_int = bg_prot_dsX;
            bg_prot_fg_er_int(fg_er_dsX==0) = 0;

            feats = zeros(1,4);
            feats(1) = bg_int / (prot_int + bg_int);
            feats(2) = sum(bg_prot_fg_nuc_int(:)) / sum(fg_nuc_dsX(:));
            feats(3) = sum(bg_prot_fg_tub_int(:)) / sum(fg_tub_dsX(:));
            feats(4) = sum(bg_prot_fg_er_int(:)) / sum(fg_er_dsX(:));

            names = {'non_object_fluorescence', ...
                'nuc:non_object_fluorescence_int_to_nuc_fg_int', 'tub:non_object_fluorescence_int_to_tub_fg_int', 'er:non_object_fluorescence_int_to_er_fg_int'};
            slfnames = {'7.79','','',''};
            
          case 'objField'
            prot = loadChannel( readlist{i}, prot);
            nuc = loadChannel( readlist_nuc{i}, nuc);
            
            prot_dsX = downsampleChannel( prot, DSF, prot_dsX);
            nuc_dsX = downsampleChannel( nuc, DSF, nuc_dsX);
            
            prot_dsX_filtered = imfilter( prot_dsX, fspecial('disk',5));
            thr1 = 255*graythresh(prot_dsX_filtered(prot_dsX_filtered>min(prot_dsX_filtered(:))));
            thr2 = 255*graythresh(prot_dsX_filtered(prot_dsX_filtered>thr1));
            fg_prot_dsX_filtered = prot_dsX_filtered;
            fg_prot_dsX_filtered(prot_dsX_filtered<thr2) = 0;
            
            [fg_prot_dsX thr_prot_dsX] = thresholdChannel( prot_dsX, fg_prot_dsX, thr_prot_dsX);
            [fg_nuc_dsX thr_nuc_dsX] = thresholdChannel( nuc_dsX, fg_nuc_dsX, thr_nuc_dsX);

            bwl_prot = bwlabel(fg_prot_dsX>0,4);
            bwl_nuc = bwlabel(fg_nuc_dsX>0,4);
            bwl_prot_filtered = bwlabel(fg_prot_dsX_filtered>0,4);

            props_nuc = regionprops(bwl_nuc,'Area');
            stats_nuc = zeros(length(props_nuc),1);
            for j = 1:length(props_nuc)
                stats_nuc(j,:) = props_nuc(j).Area;
            end
            idx_nuc = find(stats_nuc>50);
            
            props_prot = regionprops(bwl_prot,'Area');
            stats_prot = zeros(length(props_prot),1);
            for j = 1:length(props_prot)
                stats_prot(j,:) = props_prot(j).Area;
            end
            idx_prot = find(stats_prot>2);
            
            props_prot_filtered = regionprops(bwl_prot_filtered,'Area');
            stats_prot_filtered = zeros(length(props_prot_filtered),1);
            for j = 1:length(props_prot_filtered)
                stats_prot_filtered(j,:) = props_prot_filtered(j).Area;
            end
            idx_prot_filtered = find(stats_prot_filtered>2);
            
            feats = zeros(1,5);
            feats(1) = mean(stats_prot(idx_prot));
            feats(2) = var(stats_prot(idx_prot));
            feats(3) = max(stats_prot(idx_prot))/min(stats_prot(idx_prot));
            feats(4) = length(idx_prot)/length(idx_nuc);
            feats(5) = length(idx_prot_filtered)/length(idx_nuc);
            
            names = {'Avg. # pixels per obj.','Var. of # pixels per obj.','Ratio of size of largest obj. to smallest', ...
                'nuc: # prot objs to nuc objs', 'nuc: # cent objs to nuc objs'};
            slfnames = {'SLF1.3','SLF1.4','SLF1.5','',''};
                        
          case 'edgeField'
            prot = loadChannel( readlist{i}, prot);
            prot_dsX = downsampleChannel( prot, 2, prot_dsX);
            [fg_prot_dsX thr_prot_dsX] = thresholdChannel( prot_dsX, fg_prot_dsX, thr_prot_dsX);

            [names,feats,slfnames] = ml_imgedgefeatures( double(fg_prot_dsX));
            
          case 'skelField'
            prot = loadChannel( readlist{i}, prot);
            prot_dsX = downsampleChannel( prot, 2, prot_dsX);
            [fg_prot_dsX thr_prot_dsX] = thresholdChannel( prot_dsX, fg_prot_dsX, thr_prot_dsX);

            imageproc = double(fg_prot_dsX);
            imageproc = imageproc / max(imageproc(:));
            [names,feats,slfnames] = ml_imgskelfeats( imageproc);
            
          otherwise
            error('improper feature set specified');
        end
        impath = readlist{i};

        % saving results
        save( writepath,'feats','names','slfnames','impath');

        fclose(fid);
        delete(tmpfile);
    end
end

return

