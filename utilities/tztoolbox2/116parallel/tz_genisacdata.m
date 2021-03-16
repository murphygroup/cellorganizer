function tz_genisacdata(infofile,resultdir)
%TZ_GENISACDATA Generate data for ISAC.

if license('checkout','statistics_toolbox')==0
  return
end

rand('state',sum(100*clock));

addpath /home/tingz/matlab/shared
tz_initpath

load(infofile);

nCluster  = max(heladata.clusterLabels);

for i=1:heladata.nclass
    clusterLabels{i} = unique(heladata.clusterLabels( ...
        heladata.combobjclass==i))';
end

for i=1:nCluster
    clusterObjects{i} = tz_gethelainfo(heladata, ...
        struct('infotype','objects','cluster',i));
end

imgsize = [380 512];
hulls = [1 1 1 1 1 1 1 1 1 1];

for i=1:length(heladata.osm)
    for j=1:20
        
        imgfile = [heladata.classes{i} num2str(j) '.tif'];
        controlDirectory = [resultdir filesep imgfile '.ctr'];
        [s,msg] = mkdir(resultdir,[imgfile '.ctr']);
        
        %If the job has been taken over
        if strfind(msg,'exist')
          continue
        end
        resultFile = [resultdir filesep imgfile]

        if ~exist(resultFile,'file')
            [objnum,dists] = tz_genosmparam(heladata.osm(i));

            cellobjs = {};
            for k=1:length(objnum)
                if objnum(k)>0
                    genobjs = tz_genobj(clusterObjects{clusterLabels{i}(k)}, ...
                        objnum(k),struct('model','aut','hull', ...
                        hulls(i),'sizeThreshold',4000,'hullshape',0));

                    cellobjs = {cellobjs{:},genobjs{:}};
                end
            end
            img = tz_mobj2img(cellobjs,struct('method','syn','imgsize', ...
                imgsize,'dists',dists));

            imwrite(mat2gray(img),resultFile);
        end
        rmdir(controlDirectory)
    end
    
end
