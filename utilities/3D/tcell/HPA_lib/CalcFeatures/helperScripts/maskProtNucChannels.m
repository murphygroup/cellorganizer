
tmp_indexer = zeros(2,1);
if isempty( protstruct.channel_regions)
    tmp_indexer(1) = 1;
end
if isempty( nucstruct.channel_regions)
    tmp_indexer(2) = 1;
end

if sum(tmp_indexer==0)~=length(tmp_indexer)
    bwl = bwlabel(maskstruct.channel==max(maskstruct.channel(:)),4);
    s = size(maskstruct.channel);
    ma = max(bwl(:));

    for i=1:ma
        [r c] = find(bwl==i);
        minr = min(r);  maxr = max(r);  minc = min(c);  maxc = max(c);
        newIm = uint8(zeros(maxr-minr+1,maxc-minc+1));
        idx = (c-1)*s(1)+r;
        idx2 = (c-min(c))*(max(r)-minr+1)+r-min(r)+1;

        if tmp_indexer(1)
            newIm(idx2) = protstruct.channel(idx);
            protstruct.channel_regions{i} = newIm;
        end
    
        if tmp_indexer(2)
            newIm(idx2) = nucstruct.channel(idx);
            nucstruct.channel_regions{i} = newIm;
        end
    end
end
