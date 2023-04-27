function jaccard_index = align_calc_jaccard(img1,img2)
% created 3/23/2023 R.F.Murphy
% 4/7/2023 R.F. Murphy trap and report erros
% 4/13/2023 R.F. Murphy fix jaccard index calculation

try
% requires that the two images are the same size
    if(size(img1)~=size(img2))
        warning('In align_calc_jaccard: mismatched image sizes');
    else

    s1 = regionprops(img1,'centroid');
    centroid1=s1.Centroid;
    s2 = regionprops(img2,'centroid');
    centroid2=s2.Centroid;
    shft = centroid2-centroid1;
    if any(shft)
        shft = round(shft);
        img1shifted = imtranslate(img1,shft);
        centroid1shifted = regionprops(img1shifted,'centroid');
    else
        img1shifted = img1;
    end
    orimg=or(img1shifted,img2);
    andimg=and(img1shifted,img2);
    jaccard_index=sum(andimg(:))/sum(orimg(:)); %4/13/2023
    end
catch
    warning('In align_calc_jaccard: calculation failed.');
    jaccard_index=NaN;
end