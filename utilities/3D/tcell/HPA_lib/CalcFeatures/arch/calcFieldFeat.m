% function fieldFeatCalc

addpath ./lib/system/
addpath ./lib/features/
addpath ./lib/segmentation/output/
addpath /home/jnewberg/CVS/SLIC
setcvspath('/home/jnewberg/CVS/SLIC/')

fsetnames = {...
    'overlapField', ...
    'nonObjFluorField', ...
    'objField', ...
    'mutualInfo16Field', ...
    'texture64Field', ...
    'texture64x2Field', ...
    'texture64x4Field', ...
    'tasField', ...
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
    'texture256x4Field', ...
    'mutualInfo32Field', ...
    'mutualInfo64Field', ...
    'edgeField' ...
};

%     'skelField'... takes toooooo long
% fsetnames = {'overlapField'};


fsetnames = {...
    'overlapField', ...
    'nonObjFluorField', ...
    'objField', ...
    'mutualInfo16Field', ...
    'texture64Field', ...
    'texture64x2Field', ...
    'texture64x4Field', ...
    'tasField' ...
};



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


for i=1:length(readlist)
    prot = [];
    nuc = [];
    tub = [];
    er = [];
    prot_ds2 = [];
    nuc_ds2 = [];
    tub_ds2 = [];
    er_ds2 = [];
    thr_prot_ds2 = [];
    protimg_prot_ds2 = [];
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

        if isempty(prot)
            prot = imread(readlist{i});
            sprot = size(prot);
        end


        % load additional channels
        switch fsetnames{zed}
          case {'mutualInfo16Field','mutualInfo32Field','mutualInfo64Field',...
              'overlapField'}
            if isempty(nuc)
                nuc = imread(readlist_nuc{i});
            end
            if isempty(tub)
                nuc = imread(readlist_tub{i});
            end
            if isempty(er)
                nuc = imread(readlist_er{i});
            end
          otherwise
            error('Problem loading reference channels');
        end

        % downsample images
        switch fsetnames{zed}
          case {'texture16x2Field','texture32x2Field','texture64x2Field', ...
                'texture128x2Field','texture256x2Field', ...
                'overlapField', 'edgeField','nonObjectFluorField',...
                'objField', 'skelField'}
            if isempty(prot_ds2)
                prot_ds2 = imresize(prot,[sprot(1) sprot(2)]/2);
            end
            if isempty(nuc_ds2)
                if ~isempty(nuc)
                    nuc_ds2 = imresize(nuc,[sprot(1) sprot(2)]/2);
                end
            end
            if isempty(tub_ds2)
                if ~isempty(tub)
                    tub_ds2 = imresize(tub,[sprot(1) sprot(2)]/2);
                end
            end
            if isempty(er_ds2)
                if ~isempty(er)
                    er_ds2 = imresize(er,[sprot(1) sprot(2)]/2);
                end
            end
          otherwise
            error('Problem downsampling image(s)');
        end

        % threshold images
        switch fsetnames{zed}
          case {'texture16x2Field','texture32x2Field','texture64x2Field', ...
                'texture128x2Field','texture256x2Field', ...
                'overlapField', 'edgeField','nonObjectFluorField',...
                'objField', 'skelField'}
            if isempty(thr_prot_ds2)
                if ~isempty(prot_ds2)
                    thr_prot_ds2 = graythresh(prot_ds2(prot_ds2>0))*255;
                    fg_prot_ds2 = prot_ds2;
                    fg_prot_ds2(prot_ds2<=thr_prot_ds2) = 0;
                end
            end
          otherwise
            error('Problem thresholding image');
        end



        switch fsetnames{zed}
          case 'tasField'
            [names, feats, slfnames] = ml_tas( prot,0);
            
          case {'mutualInfo16Field','mutualInfo32Field','mutualInfo64Field'}
            NBINS = str2double( fsetnames{zed}(end-6:end-5));
            NBINS(2) = NBINS(1);

            prot_re = reshape(prot,[s(1)*s(2) 1]);
            nuc_re = reshape(nuc,[s(1)*s(2) 1]);
            tub_re = reshape(tub,[s(1)*s(2) 1]);
            er_re = reshape(er,[s(1)*s(2) 1]);

            feats1 = mutualinformation( [prot_re nuc_re], NBINS);
            feats2 = mutualinformation( [prot_re tub_re], NBINS);
            feats3 = mutualinformation( [prot_re er_re], NBINS);

            feats1(1,2) = corr2( prot_re, nuc_re);
            feats2(1,2) = corr2( prot_re, tub_re);
            feats3(1,2) = corr2( prot_re, er_re);

            feats = [feats1 feats2 feats3];
            str = [ '_' num2str(NBINS(1)) '_graylevels'];
            names = {['nuc:mutual_information' str],['nuc:correlation' str], ...
                ['tub:mutual_information' str],['tub:correlation' str], ...
                ['er:mutual_information' str],['er:correlation' str]};
            slfnames = repmat({''},[1 length(names)]);
            
          case {'texture16Field','texture32Field','texture64Field', ...
                'texture128Field','texture256Field'}
            GLEVELS = str2double( fsetnames{zed}(end-7:end-5));
            if isnan(GLEVELS)
                GLEVELS = str2double( fsetnames{zed}(end-6:end-5));
            end
            [feats,names,slfnames] = features_texture( prot, GLEVELS);
            
          case {'texture16x2Field','texture32x2Field','texture64x2Field', ...
                'texture128x2Field','texture256x2Field'}
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
            GLEVELS = str2double( fsetnames{zed}(end-9:end-7));
            if isnan(GLEVELS)
                GLEVELS = str2double( fsetnames{zed}(end-8:end-7));
            end
            prot_ds4 = imresize(prot, [sprot(1) sprot(2)]);
            [feats,names,slfnames] = features_texture( prot_ds4, GLEVELS);
            
            for m=1:length(names)
                names{m} = [names{m} '_4x_downsampled'];
            end
            
          case 'overlapField'
            [feats1,names1] = features_overlap( prot_ds2, nuc, 'nuc');
            [feats2,names2] = features_overlap( prot_ds2, tub, 'tub');
            [feats3,names3] = features_overlap( prot_ds2, er, 'er');

            feats = [feats1 feats2 feats3];
            names = [names1 names2 names3];
            slfnames = repmat({''},[1 length(names)]);
            
          case 'edgeField'
            [names,feats,slfnames] = ml_imgedgefeatures( double(fg_prot_ds2);
            
          case 'nonObjFluorField'
            prot_int = sum(double(fg_prot_ds2(:)));
            bg_prot_ds2 = prot_ds2;
            bg_prot_ds2(prot_ds2>thr_prot_ds2) = 0;
            bg_int = sum(double(bg_prot_ds2(:)));
            feats = bg_int / (prot_int + bg_int);
            names = 'non_object_fluorescence';
            slfnames = {'7.79'};
            
          case 'objField'
            bwl = bwlabel(fg_prot_ds2>0,4);
            props = regionprops(bwl,'Area');
            stats = zeros(length(props),1);
            for j = 1:length(props)
                stats(j,:) = props(j).Area;
            end
            stats(stats<3) = [];
            feats(1) = mean(stats);
            feats(2) = var(stats);
            feats(3) = max(stats)/min(stats);
            names = {'Avg. # pixels per obj.','Var. of # pixels per obj.','Ratio of size of largest obj. to smallest'};
            slfnames = {'SLF1.3','SLF1.4','SLF1.5'};
                        
          case 'skelField'
            imageproc = double(fg_prot_ds2);
            imageproc = imageproc / max(imageproc(:));
            [names,feats,slfnames] = ml_imgskelfeats( imageproc);
            
          otherwise
            error('improper feature set specified');
        end
        impath = readlist{i};

        % saving results
        save( writelist{i},'feats','names','slfnames','impath');





        fclose(fid);
        delete(tmpfile);
    end
end



















break
return

for zed=1:length(fsetnames)
    mout = findreplacestring( uout, '/', '_');
    mout = findreplacestring( mout, '.tif', ['_' fsetnames{zed} '.mat']);
    mout = findreplacestring( mout, rootdir_, ['./data/features/field/']);

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

        prot = imread(readlist{i});
        sprot = size(prot);

        switch fsetnames{zed}
          case 'tasField'
            [names, feats, slfnames] = ml_tas( prot,0);
            
          case {'mutualInfo16Field','mutualInfo32Field','mutualInfo64Field'}
            NBINS = str2double( fsetnames{zed}(end-6:end-5));
            NBINS(2) = NBINS(1);

            nuc = imread(readlist_nuc{i});
            tub = imread(readlist_tub{i});
            er = imread(readlist_er{i});

            prot = reshape(prot,[s(1)*s(2) 1]);
            nuc = reshape(nuc,[s(1)*s(2) 1]);
            tub = reshape(tub,[s(1)*s(2) 1]);
            er = reshape(er,[s(1)*s(2) 1]);

            feats1 = mutualinformation( [prot nuc], NBINS);
            feats2 = mutualinformation( [prot tub], NBINS);
            feats3 = mutualinformation( [prot er], NBINS);

            feats1(1,2) = corr2( prot, nuc);
            feats2(1,2) = corr2( prot, tub);
            feats3(1,2) = corr2( prot, er);

            feats = [feats1 feats2 feats3];
            str = [ '_' num2str(NBINS(1)) '_graylevels'];
            names = {['nuc:mutual_information' str],['nuc:correlation' str], ...
                ['tub:mutual_information' str],['tub:correlation' str], ...
                ['er:mutual_information' str],['er:correlation' str]};
            slfnames = repmat({''},[1 length(names)]);
            
          case {'texture16Field','texture32Field','texture64Field', ...
                'texture128Field','texture256Field'}
            GLEVELS = str2double( fsetnames{zed}(end-7:end-5));
            if isnan(GLEVELS)
                GLEVELS = str2double( fsetnames{zed}(end-6:end-5));
            end
            [feats,names,slfnames] = features_texture( prot, GLEVELS);
            
          case {'texture16x2Field','texture32x2Field','texture64x2Field', ...
                'texture128x2Field','texture256x2Field'}
            GLEVELS = str2double( fsetnames{zed}(end-9:end-7));
            if isnan(GLEVELS)
                GLEVELS = str2double( fsetnames{zed}(end-8:end-7));
            end
            s = sprot / 2;
            prot = imresize(prot, [s(1) s(2)]);
            [feats,names,slfnames] = features_texture( prot, GLEVELS);
            
            for m=1:length(names)
                names{m} = [names{m} '_2x_downsampled'];
            end

          case {'texture16x4Field','texture32x4Field','texture64x4Field', ...
                'texture128x4Field','texture256x4Field'}
            GLEVELS = str2double( fsetnames{zed}(end-9:end-7));
            if isnan(GLEVELS)
                GLEVELS = str2double( fsetnames{zed}(end-8:end-7));
            end
            s = sprot / 4;
            prot = imresize(prot, [s(1) s(2)]);
            [feats,names,slfnames] = features_texture( prot, GLEVELS);
            
            for m=1:length(names)
                names{m} = [names{m} '_4x_downsampled'];
            end
            
          case 'overlapField'
            nuc = imread(readlist_nuc{i});
            er = imread(readlist_er{i});
            tub = imread(readlist_tub{i});

            s = sprot;
            prot_ds = imresize(prot,[s(1) s(2)]/2);
            nuc_ds = imresize(nuc,[s(1) s(2)]/2);
            tub_ds = imresize(tub,[s(1) s(2)]/2);
            er_ds = imresize(er,[s(1) s(2)]/2);
            
            [feats1,names1] = features_overlap( prot_ds, nuc, 'nuc');
            [feats2,names2] = features_overlap( prot_ds, tub, 'tub');
            [feats3,names3] = features_overlap( prot_ds, er, 'er');

            feats = [feats1 feats2 feats3];
            names = [names1 names2 names3];
            slfnames = repmat({''},[1 length(names)]);
            
          case 'edgeField'
            s = sprot;
            prot_ds = imresize(prot,[s(1) s(2)]/2);

            imageproc = double(prot - thr_prot_ds2);
            [names,feats,slfnames] = ml_imgedgefeatures( imageproc);
                        
          case 'nonObjFluorField'
            s = sprot;
            prot_ds = imresize(prot,[s(1) s(2)]/2);

            prot_obj = sum(prot(:) > thr_prot_ds2);
            prot_nobj = sum(prot(:) <= thr_prot_ds2);
            feats = prot_nobj / (prot_nobj + prot_obj);
            names = 'non_object_fluorescence';
            slfnames = {'7.79'};
            
          case 'objField'
            s = sprot;
            prot_ds = imresize(prot,[s(1) s(2)]/2);

            bwl = bwlabel(prot_ds(:) > thr_prot_ds2,4);
            props = regionprops(bwl,'Area');
            stats = zeros(length(props),1);
            for j = 1:length(props)
                stats(j,:) = props(j).Area;
            end
            stats(stats<3) = [];
            feats(1) = mean(stats);
            feats(2) = var(stats);
            feats(3) = max(stats)/min(stats);
            names = {'Avg. # pixels per obj.','Var. of # pixels per obj.','Ratio of size of largest obj. to smallest'};
            slfnames = {'SLF1.3','SLF1.4','SLF1.5'};
                        
          case 'skelField'
            s = sprot;
            prot_ds = imresize(prot,[s(1) s(2)]/2);

            imageproc = double(prot_ds - thr_prot_ds2);
            imageproc = imageproc / max(imageproc(:));
            [names,feats,slfnames] = ml_imgskelfeats( imageproc);
            
          otherwise
            error('improper feature set specified');
        end
        impath = readlist{i};

        % saving results
        save( writelist{i},'feats','names','slfnames','impath');

        fclose(fid);
        delete(tmpfile);
    end
end
