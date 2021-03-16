function tz_kmeansaic_p(k)

load /home/tingz/3dobjects/combfeats.mat

[aics,centers,posts]=tz_kmeansaic(combfeats,k,10,1,0,1);

save(['/home/tingz/3dobjects/clustering1/','kmeansaic',num2str(k),'.mat'],'aics','centers','posts')
