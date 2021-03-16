function tz_findimgobjs_3d_mcf(dirname,savedir,loadthre)
%TZ_FINDIMGOBJS_3D_MCF Obsolete. See ML_FINDIMGOBJS_3D_MCF.

%function tz_findimgobjs_3d_mcf(dirname,savedir,threshold)
%   
%OVERVIEW:
%   find objects in mcf 3d image files
%PARAMETERS:
%   dirname - root image directory
%   savedir -  save directory
%   loadthre - loading threshold
%RETURN:
% 
%DESCRIPTION:
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   04-NOV-2004 Modified TINGZ
%       - add comments

error(tz_genmsg('of','tz_findimgobjs_3d_mcf','ml_findimgobjs_3d_mcf'));

class_names=tz_cleandirs(mv_dir(dirname));
% rmi=[];
% for i=1:length(class_names)
%     if class_names{i}(1)=='.'
%         rmi=[rmi,i];
%     end
% end
% 
% class_names(rmi)=[];

if ~exist(savedir,'dir')
    unix(['mkdir ' savedir]);
end


NumberOfClasses = length(class_names);
first_class = 1;
last_class = NumberOfClasses;
for class = first_class : last_class
    class_name = class_names{class};
    status=mkdir(savedir,class_name);
    saveclassdir=[savedir '/' class_name];
    
    cells=tz_cleandirs(mv_dir([dirname '/' class_name]));
    ncell=length(cells);
    for N=1:ncell
       % status=mkdir(saveclassdir,cells{N});
        %savecelldir=[saveclassdir '/' cells(N)];
        celldir=[dirname '/' class_name '/' cells{N}];
        protdir = [celldir '/prot'];
        cropfiles = mv_dir([celldir '/crop/*.tif*']);
        objects=tz_findimgobjs_3d(protdir,[celldir '/crop/' cropfiles{1}],loadthre);
        save([saveclassdir '/' cells{N} '.mat'],'objects')
    end
end

function cnames=tz_cleandirs(names)

rmi=[];
for i=1:length(names)
    if names{i}(1)=='.'
        rmi=[rmi,i];
    end
end

cnames=names;
cnames(rmi)=[];
