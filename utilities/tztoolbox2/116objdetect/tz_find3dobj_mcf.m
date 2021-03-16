function tz_find3dobj_mcf(rootdir,savedir,channel)
%TZ_FIND3DOBJ_MCF Obsolete.

%function tz_find3dobj_mcf(rootdir,savedir)

error(tz_genmsg('of','tz_find3dobj_mcf'));

if ~exist('channel','var')
    channel='prot';
end

classes=tz_cleandirs(mv_dir(rootdir));

nclass=length(classes);


for i=1:nclass
    protdir=[rootdir '/' classes{i}];
    cells=tz_cleandirs(mv_dir(protdir));
    ncell=length(cells);
    status=mkdir(savedir,classes{i});
    
    for j=1:ncell
        celldir=[protdir '/' cells{j}];
        
        cropdir=[celldir '/' 'crop'];
        crops=tz_cleandirs(mv_dir(cropdir));
        ncrop=length(crops);
        imgloaded=0;
        [i j]
        for k=1:ncrop
            savefile=[cells{j} 'crop' num2str(k) '.mat'];
            saveprotdir=[savedir '/' classes{i}];
            if ~exist([saveprotdir '/' savefile],'file')
                if imgloaded==0
                    protimg=tz_loadimage([celldir '/' channel],'tif',[]);
                    imgloaded=1;
                end
                cropimg=imread([cropdir '/' crops{k}]);
                
                objects=tz_find3dobj(protimg,cropimg,1);
                save([saveprotdir '/' savefile],'objects');
            end
        end
    end
end