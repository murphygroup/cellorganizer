function learngaussobjs_run(imgdir,savedir,param)
% Fit gaussian distributions to the objects

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param, ...
    struct('modelname','vesicle','imageSize',[],'disp',0));
model.name = param.modelname;

if ~isfield(param,'ml_objgaussmix')
    param.ml_objgaussmix = struct([]);
end

param.ml_objgaussmix = ml_initparam(param.ml_objgaussmix, ...
    struct('filter',fspecial3('Gaussian',[10,10,3]),'mindist',5, ...
    'isshow',0,'gmm',struct('covartype','full',...
    'options',[0 1 0 0 0 0 0 0 0 0 0 0 0 1000 1 20])));

for i = 1:50
    load(['./inter_results/protein_objects_gaussian/original_objects/' protype '/obj' num2str(i) '.mat'])
    disp([protype ' image ' num2str(i)])   
    %learn Gaussian mixture model
    ct = 0;
    for r = 1:10
        for c = 1:10
            blockobjs = objects{r,c};
            if ~isempty(blockobjs)
                disp(['Block (' int2str(r) ', ' int2str(c) ')'])
                for j = 1:length(blockobjs)
                    ct = ct + 1;
                    obj = blockobjs{j};

                    %ERROR
                    objintens{i}(ct) = sum(obj(:,end));
                    mixes{i}(ct) = ml_objgaussmix(obj,...
                        [],param.ml_objgaussmix);
                    offsets{i}(ct,:) = 100*[r-1,c-1,0];


                end
            end
        end
    end
    save([savedir '/' protype '_gaussobjs'],...
        'mixes','objintens','offsets')
end
