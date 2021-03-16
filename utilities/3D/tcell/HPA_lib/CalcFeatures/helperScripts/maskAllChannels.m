
if isempty(protfieldstruct.channel)
    protfieldstruct.channel = imread(protfieldstruct.channel_path);
end

if isempty(nucfieldstruct.channel)
    nucfieldstruct.channel = imread(nucfieldstruct.channel_path);
end

if isempty(tubfieldstruct.channel)
    tubfieldstruct.channel = imread(tubfieldstruct.channel_path);
end

if isempty(erfieldstruct.channel)
    erfieldstruct.channel = imread(erfieldstruct.channel_path);
end

if isempty(maskfieldstruct.channel)
    maskfieldstruct.channel = imread(maskfieldstruct.channel_path);
end


tmp_indexer = zeros(4,1);
if isempty( protfieldstruct.channel_regions)
    tmp_indexer(1) = 1;
end
if isempty( nucfieldstruct.channel_regions)
    tmp_indexer(2) = 1;
end
if isempty( tubfieldstruct.channel_regions)
    tmp_indexer(3) = 1;
end
if isempty( erfieldstruct.channel_regions)
    tmp_indexer(4) = 1;
end

if sum(tmp_indexer==0)~=length(tmp_indexer)
    bwl = bwlabel(maskfieldstruct.channel==max(maskfieldstruct.channel(:)),4);
    s = size(maskfieldstruct.channel);
    ma = max(bwl(:));

    for tmp_i=1:ma
        [r c] = find(bwl==tmp_i);
        minr = min(r);  maxr = max(r);  minc = min(c);  maxc = max(c);
        newIm = uint8(zeros(maxr-minr+1,maxc-minc+1));
        idx = (c-1)*s(1)+r;
        idx2 = (c-min(c))*(max(r)-minr+1)+r-min(r)+1;

        if tmp_indexer(1)
            newIm(idx2) = protfieldstruct.channel(idx);
            protfieldstruct.channel_regions{tmp_i} = newIm;
        end
    
        if tmp_indexer(2)
            newIm(idx2) = nucfieldstruct.channel(idx);
            nucfieldstruct.channel_regions{tmp_i} = newIm;
        end

        if tmp_indexer(3)
            newIm(idx2) = tubfieldstruct.channel(idx);
            tubfieldstruct.channel_regions{tmp_i} = newIm;
        end
    
        if tmp_indexer(4)
            newIm(idx2) = erfieldstruct.channel(idx);
            erfieldstruct.channel_regions{tmp_i} = newIm;
        end
    end
end
