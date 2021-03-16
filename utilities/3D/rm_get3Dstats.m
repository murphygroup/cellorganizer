function stats = rm_get3Dstats(img)  

stats.Volume = 0;
stats.SurfaceArea = 0;
stats.Eccentricity = NaN;
stats.MajorAxisLength = 0;
MaxArea = 0;
for k=1:size(img,3)
    nucbin = (img(:,:,k)==255);
    if sum(nucbin(:))>0
        s=regionprops(nucbin,{'Area', 'Perimeter', 'Eccentricity', 'MajorAxisLength'});
        stats.Volume = stats.Volume + s.Area;
        stats.SurfaceArea = stats.SurfaceArea + s.Perimeter;
        % save Eccentricity and MajorAxisLength only for the biggest slice
        if s.Area > MaxArea
            MaxArea = s.Area;
            stats.Eccentricity = s.Eccentricity;
            stats.MajorAxisLength = s.MajorAxisLength;
        end
    end
end


