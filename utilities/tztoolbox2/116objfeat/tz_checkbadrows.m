function [ncellobjc,ncellobjf]=tz_checkbadrows(all_coords,all_features)
%TZ_CHECKBADROWS Check if two sets match with each other on object numbers.
%   NCELLOBJC = TZ_CHECKBADROWS(ALL_COORDS,ALL_FEATURES)
%   
%   [NCELLOBJC,NCELLOBJF] = TZ_CHECKBADROWS(...)
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function [ncellobjc,ncellobjf]=tz_checkbadrows(all_coords,all_features)
%
%OVERVIEW:
%   check if the object numbers in two data match or not
%PARAMETERS:
%   all_coords - mcf objects
%   all_features - mcf object features
%RETURN:
%   ncellobjc - object number in each cell (all_coords)
%   ncellobjf - object number in each cell (all_features)
%HISTORY:
%   ??-???-???? Initial write TINGZ
%   03-NOV-2004 Modified TINGZ
%       - add comments

nclass=length(all_coords);
nobjc=zeros(1,10);
nobjf=zeros(1,10);
ncellobjc=zeros(10,100);
ncellobjf=zeros(10,100);

for i=1:nclass
    ncellc(i)=length(all_coords{i});
    ncellf(i)=length(all_features{i});
    
    for j=1:ncellc(i)
        nobj=size(all_coords{i}{j},1);
        ncellobjc(i,j)=nobj;
        nobjc(i)=nobjc(i)+nobj;
    end       
    for j=1:ncellf(i)
        nobj=size(all_features{i}{j},1);
        ncellobjf(i,j)=nobj;
        nobjf(i)=nobjf(i)+nobj;
    end
end

[ncellc;ncellf;ncellc-ncellf;nobjc;nobjf;nobjc-nobjf]