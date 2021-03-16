function tz_kmeansaic2_p(k)

tz_pinitpath
load pp_combobjfeats.mat
rand('state',[1:35]'+2);
[aics,centers,posts]=tz_kmeansaic(combobj(:,[1 3 5:end]),k,10,1,0,0);

save(['/home/tingz/paper2/kvsaic/','kmeansaic',num2str(k),'.mat'], ...
    'aics','centers','posts')
