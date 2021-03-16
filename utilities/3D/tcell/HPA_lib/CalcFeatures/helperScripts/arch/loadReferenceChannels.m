function [nuc,tub,er] = loadReferenceChannels( nucpath, tubpath, erpath, ...
                               nuc, tub, er)

if isempty(nuc)
    nuc = imread(nucpath);
end

if isempty(tub)
    tub = imread(tubpath);
end

if isempty(er)
    er = imread(erpath);
end

return