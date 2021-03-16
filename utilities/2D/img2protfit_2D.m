function protfit = img2protfit_2D(procimage, options_gaussmix)

if ~exist('options_gaussmix', 'var')
    options_gaussmix = []; 
end


options_gaussmix = ml_initparam(options_gaussmix, ...
                        struct('filter',fspecial('gaussian',10,10), ...
                                'mindist',10, ...
                                'isshow',0, ...
                                'gmm',struct('covartype','spherical',...
                                             'options',[0 1 0 0 0 0 0 0 0 0 0 0 0 500] ...
                                             ) ...
                                ));


latents = [];
intensities = [];

disp('Learning Gaussian mixtures');
tic
% load(procfile)
objs = ml_findobjs2D(procimage);
toc

objcof = [];
for j=1:length(objs)
    disp(['Object:' num2str(j)]);
    obj = objs{j};
    objintens(j) = sum(obj(:,3));
    if size(obj,1)==1
        mixes{j} = [];
    else
        tic
        mixture = ml_objgaussmix2D( obj,options_gaussmix);
        mixes{j} = mixture;

        if options_gaussmix.isshow == 1
            drawnow
        end

        toc

        %Get positions

        objintensity = objintens(j);
        if ~isempty(mixes{j})
            for k=1:mixes{j}.ncentres
                latents = [latents; mixes{j}.covars(1,1,k)];
                intensities = [intensities; ...
                               objintensity*mixes{j}.priors(k)];
                objcof(end+1,:) = mixes{j}.centres(k,:);
            end
        end
   
    end
end

%remove illegal entries
rminds = isnan(latents) | isnan(intensities);
latents(rminds) = [];
intensities(rminds) = [];

protfit.objs = objs;
protfit.mixes = mixes;
protfit.latents = latents;
protfit.intensities = intensities;
protfit.objcof = objcof;
protfit.objnums = size(objcof,1);
protfit.options = options_gaussmix;
