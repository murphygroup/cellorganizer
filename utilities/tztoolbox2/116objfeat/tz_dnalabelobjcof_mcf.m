function lobjcofs=tz_dnalabelobjcof_mcf(objcofs,rootdir)
%TZ_DNALABELOBJCOF_MCF Label object COF according to nucleus overlap.
%   LOBJCOFS = TZ_DNALABELOBJCOF_MCF(OBJCOFS,ROOTDIR)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function lobjcofs=tz_dnalabelobjcof_mcf(objcofs,rootdir)
%OVERVIEW:
%   label object COF by 1 if it is in the nucleus or 0 if not
%PARAMETERS:
%   objcofs - object COF
%   rootdir - directory storing nuc edges
%RETURN
%   lobjcofs - labeled COF
%DESCRIPTION
%
%HISTORY:
%   27-SEP-2004 Initial write TINGZ

class_names=tz_cleandirs(mv_dir(rootdir));

NumberOfClasses = length(class_names);
first_class = 1;
last_class = NumberOfClasses;

lobjcofs=objcofs

for class = first_class : last_class
    class_name = class_names{class};

    dnadir=[rootdir '/' class_name '/dnaedge'];
    
    dnas=tz_cleandirs(mv_dir(dnadir));

    ndna=length(dnas);
    for N=1:ndna
        load([dnadir '/' dnas{N}]);
        ledge=bwlabel(dnaedge);
        objnum=max(ledge(:));
        
        lhist=[];
        for i=1:objnum
            lhist(i)=sum(sum(ledge==i));
        end
        [y,maxl]=max(lhist);
        dnaedge(ledge~=maxl)=0;
        
        bwnuc=imfill(dnaedge,'hole');
        nobj=size(objcofs{class}{N},1);
        labels=[];
        for i=1:nobj
            cof=round(objcofs{class}{N}(i,:));
            labels(i)=bwnuc(cof(1),cof(2));
        end
        lobjcofs{class}{N}=[lobjcofs{class}{N},labels'];
    end
end