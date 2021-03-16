function objects=tz_findmainobj_bw_mcf(rootdir,srcdir)

%function objects=tz_findmainobj_bw_mcf(rootdir,srcdir)
%   
%OVERVIEW:
%   find biggest objects
%PARAMETERS:
%   rootdir - root path for the images
%   srcdir - dir under class dir
%RETURN:
%   objects - 
%DESCRIPTION:
%
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   31-OCT-2004 Modified TINGZ
%       - add comments

class_names=tz_cleandirs(mv_dir(rootdir));

NumberOfClasses = length(class_names);
first_class = 1;
last_class = NumberOfClasses;

for class = first_class : last_class
    class_name = class_names{class};

    procdir=[rootdir '/' class_name '/' srcdir];
    
    imgfiles=tz_cleandirs(mv_dir(procdir));

    nimg=length(imgfiles);
    objs={};
    for N=1:nimg
        img=load([procdir '/' imgfiles{N}]);
        field=fieldnames(img);
        img=getfield(img,field{1});
        objs{N}=tz_findmainobj_bw(img);
    end
    objects{class}=objs;
end