#!/bin/bash
alias qarray='perl /export/home/install/bin/r2/qarray.pl'

declare -i i
i=$1

while [ $i -le $2 ]
do
  echo $i
  #qarray -m "tz_kmeansaicbic_p(%%)" "'/home/tingz/matlab/juggernaut/data/yeast/cellfeatall.mat','/home/tingz/matlab/juggernaut/data/yeast/clustering/cellfeat/kmeansaic',$i"
  #qarray -m "tz_kmeansaicbic_p(%%)" "'/home/tingz/matlab/juggernaut/data/helaobjfeatsforclust.mat','/home/tingz/matlab/juggernaut/data/genmodel/2DHeLa/objcluster3',$i"
  #qarray -m "tz_kmeansaicbic2_p(%%)" "'/home/tingz/matlab/juggernaut/data/yeast/cellfeatall.mat','/home/tingz/matlab/juggernaut/data/yeast/clustparam1.mat',$i"  
  #qarray -m "tz_kmeansaicbic2_p(%%)" "'/home/tingz/matlab/juggernaut/data/yeast/fieldfeatall.mat','/home/tingz/matlab/juggernaut/data/yeast/clustparam5.mat',$i"
  #qarray -m "tz_kmeansaicbic4_p(%%)" "'/home/tingz/matlab/juggernaut/data/yeast/ucsffieldfeat.mat','/home/tingz/matlab/juggernaut/data/yeast/clustparam2.mat','/home/tingz/matlab/juggernaut/data/yeast/clustering/ucsf1',$i"
  #qarray -m "tz_kmeansaicbic4_p(%%)" "'/home/tingz/matlab/juggernaut/data/yeast/cellfeat1.mat','/home/tingz/matlab/juggernaut/data/yeast/clustparam2.mat','/home/tingz/matlab/juggernaut/data/yeast/clustering/cellfeat/kmeans3',$i"
  qarray -m "tz_kmeansaicbic5_p(%%)" "'/home/tingz/matlab/juggernaut/data/yeast/cellfeat1.mat','/home/tingz/matlab/juggernaut/data/yeast/clustparam2.mat','/home/tingz/matlab/juggernaut/data/yeast/clustering/cellfeat/kmeans3',21:100"
  #qarray -m "tz_kmeansaicbic4_p(%%)" "'/home/tingz/matlab/juggernaut/data/genmodel/2DHeLa/combobjfeat.mat','/home/tingz/matlab/juggernaut/data/genmodel/2DHeLa/clustparam1.mat','/home/tingz/matlab/juggernaut/data/genmodel/2DHeLa/objcluster1',$i"
  #qarray -m "tz_kmeansaicbic4_p(%%)" "'/home/tingz/matlab/juggernaut/data/genmodel/2DHeLa/oldcombobjfeat11.mat','/home/tingz/matlab/juggernaut/data/genmodel/2DHeLa/clustparam2.mat','/home/tingz/matlab/juggernaut/data/genmodel/2DHeLa/objcluster5',$i"
  #qarray -m "tz_kmeansaicbic5_p(%%)" "'/home/tingz/matlab/juggernaut/data/genmodel/2DHeLa/combobjfeat42.mat','/home/tingz/matlab/juggernaut/data/genmodel/2DHeLa/clustparam1.mat','/home/tingz/matlab/juggernaut/data/genmodel/2DHeLa/objcluster6',81:160"
  #qarray -m "tz_kmeansaicbic5_p(%%)" "'/home/tingz/matlab/juggernaut/data/genmodel/2DHeLa/combobjfeat42.mat','/home/tingz/matlab/juggernaut/data/genmodel/2DHeLa/clustparam2.mat','/home/tingz/matlab/juggernaut/data/genmodel/2DHeLa/objcluster7',1:80"
  #qarray -m "tz_kmeansaicbic5_p(%%)" "'/home/tingz/matlab/juggernaut/data/genmodel/2DHeLa/combobjfeat11.mat','/home/tingz/matlab/juggernaut/data/genmodel/2DHeLa/clustparam3_2.mat','/home/tingz/matlab/juggernaut/data/genmodel/2DHeLa/objcluster_seeds2',1:80"
  i=$i+$3
done

echo "done"
