function [imfunc, minwinsize, maxdims, image_output_size] = get_diffeo_image_function(param)
%grj 7/25/13 Image function for Taraz's diffeomorphic model for CO
%grj 8/14/13 ability to perform only on dna segmentations, and some
%            image registration
%grj 8/29/13 Added a temp file, so we dont have to load and rotate all the
%            images to find the large bounding box
% Aug 30, 2013 G. Johnson    Changed they way files are input into the
%                            diffeomorphic model function 
    
    % xruan 01/05/2015
    if isfield(param.model.diffeomorphic, 'regfunc')
        param.model.regfunc = param.model.diffeomorphic.regfunc;
    else
        param.model = ml_initparam(param.model, struct( ...
            'regfunc', @(x,y) cellfile2registeredimg([param.preprocessingFolder filesep 'cell' num2str(x) '.mat'], y) ...
            ));
    end
    tempdir = param.model.diffeomorphic.tempdir;
%     com_align = param.model.diffeomorphic.com_align;
%     downsample = param.model.diffeomorphic.downsample;
%     bottom_align = param.model.diffeomorphic.z_align;
        
    if ~exist(tempdir, 'dir')
	   mkdir(tempdir)
    end

    imfuncparams_file = [tempdir filesep 'imfunc_params.mat'];
    numimgs = param.documentation.numimgs;

    if ~exist(imfuncparams_file, 'file')
        
        %this loop figures how how big the maximum image size after
        %registration

        imgs = cell(1,numimgs);
        regparam = cell(1,numimgs);
        imsizes = zeros(numimgs,3);
        imcrops = zeros(numimgs,4);

        for i = 1:numimgs
            disp([num2str(i) filesep num2str(numimgs)]);
            
            try
                if ~isempty(param.model.regfunc)
                    [img, regparam{i}] = param.model.regfunc(i, param.model.diffeomorphic);
                end
            catch
                img = [];
                regparam{i} = [];
                
            end
                
            if ~isempty(img)
                try
                    [imgcompress, imsize, croprange] = diffeo_compress_img(img);
                    imgs{i} = imgcompress;
                    imsizes(i,:) = imsize;
                    imcrops(i,:) = croprange(1:4);
                    
                catch
                    disp(['WARNING: Insufficient heap space to compress image ' num2str(i) '. Consider downsampling or increasing heap space.'])
                end
            end
        end
        
        maxsize = max(imsizes,[],1);
        maxsize(1:2) = max(maxsize(1:2));

        minsize = min(imsizes(~any(imsizes == 0,2),:),[],1);
        minsize(1:2) = max(minsize(1:2));
        minwinsize = ceil(maxsize(1) / ceil(maxsize(1)/minsize(1)));
        minwinsize = [minwinsize, minwinsize, minsize(3)];
        image_output_size = maxsize;
        
        maxdims = max(imsizes,[],1);
        %compress all images and throw them into the image functon
    
        imfunc = @(x) diffeo_img_function(x, imgs, image_output_size, imsizes, imcrops);
  
        save(imfuncparams_file, 'minwinsize', 'image_output_size', 'imfunc', 'maxdims', 'regparam', 'imsizes', 'imcrops', 'numimgs', '-v7.3');
    else
        load(imfuncparams_file)
    end
end


function [img] = diff_img_function(filenum, compressed_imgs, maxsize, imsizes, imcrops)
    %%% This exists for legacy purposes only.
    if ~isempty(compressed_imgs{filenum})
        img = CompressLib.decompress(compressed_imgs{filenum});
        
        img = padarray(img, [imcrops(filenum,1)-1, imcrops(filenum,3)-1], 'pre');
        img = padarray(img, [imsizes(filenum,1) - imcrops(filenum,2), imsizes(filenum,2) - imcrops(filenum,4)], 'post');

    else
        img = [];
    end
    
    if ~isempty(img)
        imsize = size(img);
        imsize = [imsize, ones(1,3 - length(imsize))];

        pad = (maxsize - imsize) ./ 2;

        img = padarray(img, floor(pad), 0, 'both');

        img = padarray(img, double(~(floor(pad) == pad)), 0, 'pre');
        img = double(img);
    end
end


