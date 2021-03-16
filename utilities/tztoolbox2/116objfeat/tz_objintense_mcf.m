function intens=tz_objintense_mcf(objects)
%TZ_OBJINTENSE_MCF Get total intensitiy of each object in a cell array.
%   INTENS = TZ_OBJINTENSE_MCF(OBJECTS) returns a [2-level cell array] in which
%   each element is a vector of intensities of the corrpesponding objects in 
%   [3-level cell array] of [object]s.
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

nclass=length(objects);

for i=1:nclass
    ncell=length(objects{i});
    for j=1:ncell
        nobjects=length(objects{i}{j});
        for k=1:nobjects
            obj=objects{i}{j}{k};
            intens{i}{j}(k,2)=sum(obj(:,3));
            intens{i}{j}(k,1)=size(obj,1);
        end
    end
end
