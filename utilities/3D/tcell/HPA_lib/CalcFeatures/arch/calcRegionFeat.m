function regionFeatCalc

addpath ./lib/system/
addpath ./lib/features/
addpath ./lib/segmentation/output/
addpath ./lib/segmentation
addpath ./lib/segmentation/seeding
addpath /home/jnewberg/CVS/SLIC
setcvspath('/home/jnewberg/CVS/SLIC/')


fsetnames = {...
    'overlapRegion', ...
    'nonObjFluorRegion', ...
    'objRegion', ...
    'mutualInfo16Region', ...
    'texture64x2Region', ...
    'texture64x4Region', ...
    'tasRegion', ...
    'edgeRegion', ...
    'skelRegion', ...
    'texture16Region', ...
    'texture32Region', ...
    'texture128Region', ...
    'texture256Region', ...
    'texture16x2Region', ...
    'texture32x2Region', ...
    'texture128x2Region', ...
    'texture256x2Region', ...
    'texture16x4Region', ...
    'texture32x4Region', ...
    'texture64x4Region', ...
    'texture128x4Region', ...
    'texture256x4Region', ...
    'mutualInfo32Region', ...
    'mutualInfo64Region' ...
};

fsetnames = {'nuclearRegion'};


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


mout = findreplacestring( uout, '/', '_');
mout = findreplacestring( mout, '.tif', '.png');
mout = findreplacestring( mout, rootdir_, ['./data/masks/']);
readlist_mask = listmatlabformat( mout);

for zed=1:length(fsetnames)
    mout = findreplacestring( uout, '/', '_');
    mout = findreplacestring( mout, '.tif', ['_' fsetnames{zed} '.mat']);
    mout = findreplacestring( mout, rootdir_, ['./data/features/region/']);

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
        nuc = imread(readlist_nuc{i});
        tub = imread(readlist_tub{i});
        er = imread(readlist_er{i});
        masks = imread(readlist_mask{i});


        seedpath = findreplacestring( readlist_mask{i},'masks','seeds');
        if strcmp( fsetnames{zed},'nuclearRegion')
            if ~exist(seedpath,'file')
                IMAGEPIXELSIZE = 0.05; % um/px
                MINNUCLEUSDIAMETER = 6; %um
                MAXNUCLEUSDIAMETER = 14; %um
                minarea = (MINNUCLEUSDIAMETER/IMAGEPIXELSIZE/2)^2*pi;
                maxarea = (MAXNUCLEUSDIAMETER/IMAGEPIXELSIZE/2)^2*pi;

                % determine foreground seeds
                fgs = fgseeds( nuc, MINNUCLEUSDIAMETER, MAXNUCLEUSDIAMETER, IMAGEPIXELSIZE);
                imwrite(fgs, seedpath);
            else
                fgs = imread(seedpath);
            end

            [masks_crop,seeds_crop] = masks2cropimgs( masks, masks, fgs);
        end


        [prot_crop,nuc_crop] = masks2cropimgs( masks, prot, nuc);
        [tub_crop,er_crop] = masks2cropimgs( masks, tub, er);
        
        allfeats = [];
        for k=1:length(prot_crop)
          prot = prot_crop{k};
          nuc = nuc_crop{k};
          tub = tub_crop{k};
          er = er_crop{k};
          s = size(prot);


          if strcmp( fsetnames{zed},'nuclearRegion')
              fgs = seeds_crop{k};
          end
          
          
          switch fsetnames{zed}
            case 'tasRegion'
              [names, feats, slfnames] = ml_tas( prot,0);
            case {'mutualInfo16Region','mutualInfo32Region','mutualInfo64Region'}
              NBINS = str2double( fsetnames{zed}(end-7:end-6));
              NBINS(2) = NBINS(1);

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
              
            case {'texture16Region','texture32Region','texture64Region', ...
                  'texture128Region','texture256Region'}
              GLEVELS = str2double( fsetnames{zed}(end-8:end-6));
              if isnan(GLEVELS)
                  GLEVELS = str2double( fsetnames{zed}(end-7:end-6));
              end
              [feats,names,slfnames] = features_texture( prot, GLEVELS);
              
            case {'texture16x2Region','texture32x2Region','texture64x2Region', ...
                  'texture128x2Region','texture256x2Region'}
              GLEVELS = str2double( fsetnames{zed}(end-10:end-8));
              if isnan(GLEVELS)
                  GLEVELS = str2double( fsetnames{zed}(end-9:end-8));
              end
              s = size(prot) / 2;
              prot = imresize(prot, [s(1) s(2)]);
              [feats,names,slfnames] = features_texture( prot, GLEVELS);
              
              for m=1:length(names)
                  names{m} = [names{m} '_2x_downsampled'];
              end
              
            case {'texture16x4Region','texture32x4Region','texture64x4Region', ...
                  'texture128x4Region','texture256x4Region'}
              GLEVELS = str2double( fsetnames{zed}(end-10:end-8));
              if isnan(GLEVELS)
                  GLEVELS = str2double( fsetnames{zed}(end-9:end-8));
              end
              s = size(prot) / 4;
              prot = imresize(prot, [s(1) s(2)]);
              [feats,names,slfnames] = features_texture( prot, GLEVELS);

              for m=1:length(names)
                  names{m} = [names{m} '_4x_downsampled'];
              end

            case 'overlapRegion'
              [feats1,names1] = features_overlap( prot, nuc);
              [feats2,names2] = features_overlap( prot, tub);
              [feats3,names3] = features_overlap( prot, er);

              feats = [feats1 feats2 feats3];
              names = [names1 names2 names3];
              slfnames = repmat({''},[1 length(names)]);

            case 'edgeRegion'
              thr = graythresh(prot)*255;
              imageproc = prot - thr;
              [names,feats,slfnames] = ml_imgedgefeatures( imageproc);

            case 'nonObjFluorRegion'
              thr = graythresh(prot)*255;
              prot_obj = sum(prot(:) > thr);
              prot_nobj = sum(prot(:) <= thr);
              feats = prot_nobj / (prot_nobj + prot_obj);
              names = 'non_object_fluorescence';
              slfnames = {'7.79'};

            case 'objRegion'
              thr = graythresh(prot)*255;
              imageproc = prot>thr;

              thr = graythresh(nuc)*255;
              dnaproc = nuc>thr;
              dnaproc = medfilt2( dnaproc, [10 10]);

              [names, feats, slfnames] = ml_imgfeatures(imageproc, dnaproc)

            case 'skelRegion'
              thr = graythresh(prot)*255;
              imageproc = double(prot - thr);
              imageproc = imageproc / max(imageproc(:));
              [names,feats,slfnames] = ml_imgskelfeats( imageproc);
 
            case 'nuclearRegion'
              [c b] = imhist(nuc);
              [a ind] = max(c);
              dnaproc = nuc - b(ind);
              
              [feats,names,slfnames] = hpa_teb_nuclear_features( dnaproc, fgs);
 
            case 'spatialRegion'
              protproc = prot > 255*graythresh(prot(prot>0));
              nucproc = nuc > 255*graythresh(nuc(nuc>0));
              erproc = er > 255*graythresh(er(er>0));
              tubproc = tub > 255*graythresh(tub(tub>0));

              DSF = 8;
              
              cellproc = tubproc|nucproc|erproc;
              [p1 names1] = spatialFeatures( protproc, nucproc, cellproc, DSF);
              [p2 names2] = spatialFeatures( protproc, erproc, cellproc, DSF);
              [p3 names3] = spatialFeatures( protproc, tubproc, cellproc, DSF);

              feats = [p1 p2 p3];
              slfnames = [names1 names2 names3];
              
              names = [];

            otherwise
              error('improper feature set specified');
            end
        
          allfeats = [allfeats; feats];
        end
        
        impath = readlist{i};

        % saving results
        save( writelist{i},'allfeats','names','slfnames','impath');

        fclose(fid);
        delete(tmpfile);
    end
end
