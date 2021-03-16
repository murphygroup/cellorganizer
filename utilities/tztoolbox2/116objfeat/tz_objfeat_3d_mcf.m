function features=tz_objfeat_3d_mcf(objdir,dnadir,savedir,pclass)
%TZ_OBJFEAT_3D_MCF Under construction.

%function tz_objfeat_3d_mcf(objdir,dnadir,imgsize)

classes=tz_cleandirs(mv_dir(objdir));

nclass=length(classes);

if ~exist('pclass','var')
    pclass=1:nclass;
end

for i=pclass
    protdir=[objdir '/' classes{i}];
    protdnadir=[dnadir '/' classes{i}];
%     protimgdir=[imgdir '/' classes{i}]
    cells=tz_cleandirs(mv_dir(protdir));
    ncell=length(cells);
    status=mkdir(savedir,classes{i});
    for j=1:ncell
        [i j]
        
        protobjfile=[protdir '/' cells{j}];
        dnaobjfile=[protdnadir '/' cells{j}];
        savefile=[savedir '/' classes{i} '/' cells{j}];
        
        if ~exist(savefile,'file')
            
            protobj=load(protobjfile);
            dnaobj=load(dnaobjfile);
            
            
            for k=1:length(dnaobj.objects)
                sdna(k)=dnaobj.objects{k}.size;
            end
            [tmp,maxk]=max(sdna);
            dnaobj=[double(dnaobj.objects{maxk}.voxels);double(dnaobj.objects{maxk}.gray)]';
            
            if ~exist(savefile,'file') 
                feats=tz_objfeat_3d(protobj.objects,dnaobj,protobj.objects{1}.imgsize);
            end
            save(savefile,'feats');
        else
            load(savefile);
            
        end
        features{i}{j}=feats;
    end
end