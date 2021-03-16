function coords=tz_calcobjcof_mcf(objects)

%function coords=tz_calcobjcof_mcf(objects)
%
%OVERVIEW
%   Calculate the coordinates of COF of mcf objects
%PRAMETERS
%   objects - objects stored as mcf structure
%RETURN
%   coords - the coordinates of objects' COF
%

nclass = length(objects);

for i=1:nclass
    ncell=length(objects{i});
    cellcoords=[];
    for j=1:ncell
        cellcoords(j,:)=tz_calcobjcof(objects{i}{j});
    end
    coords{i}=cellcoords;
end
