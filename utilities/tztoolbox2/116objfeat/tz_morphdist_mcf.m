function [normdists,nucdists,celdists]=tz_morphdist_mcf(lobjcofs,rootdir,option)
%TZ_MORPHDIST_MCF Calculate normalized distances of object COFs for MCF.
%   NORMDISTS = TZ_MORPHDIST_MC(LOBJCOFS,ROOTDIR,OPTION)
%   
%   [NORMDISTS,NUCDISTS,CELLDISTS] = TZ_MORPHDIST_MC(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function [normdists,nucdists,celdists]=tz_morphdist_mcf(lobjcofs,rootdir,option)
%
%OVERVIEW
%   calculate cell morph distance
%PARAMETERS:
%   lobjcofs - labeled object cof (in nucleus 1, out 0)
%   rootdir - mcf root directory
%   option - way of normalization
%RETURN:
%   normdists - normalized distance
%   nucdists - distance to nucleus edge and cof
%   celdists - distance to cell edge
%DESCRIPTION:
%       
%HISTORY:
%   ??-OCT-2004 Initial write TINGZ
%   31-NOV-2004 Modified TINGZ


class_names=tz_cleandirs(mv_dir(rootdir));

NumberOfClasses = length(class_names);
first_class = 1;
last_class = NumberOfClasses;

normdists=lobjcofs;

for class = first_class : last_class
    class_name = class_names{class};

    dnadir=[rootdir '/' class_name '/dnaedge'];
    celldir=[rootdir '/' class_name '/celledge'];
    
    cells=tz_cleandirs(mv_dir(celldir));
    dnas=tz_cleandirs(mv_dir(dnadir));

    ndna=length(dnas);
    for N=1:ndna
        load([dnadir '/' dnas{N}]);
        load([celldir '/' cells{N}]);
        
        [normdists{class}{N},nucdists{class}{N},celdists{class}{N}]=...
            tz_morphdist(lobjcofs{class}{N},celedge,dnaedge,option);
   
    end
end