function feats=tz_loadfeats(dirname,label)
%TZ_LOADFEATS Load features from a directory.
%   FEATS = TZ_LOADFEATS(DIRNAME,LABEL) loads features from a directory
%   DIRNAME by specifying file names through LABEL, which is a structure
%   with two fields:
%       prefix: a string
%       label: a integer vector
%   The name of the ith file will be [prefix label(i)].
%   For 3D3T3 images, we have the the following naming scheme:
%       ['s' month day '-' year '-' dishlabel cell '.mat'] is
%       for static image.
%  
%       ['ts' month day '-' year '-' dishlabel cell '.mat'] is
%       for time series image.
%
%   Example:
%       if label.prefix='s1205-03-d', label.num=[1 13 114], the file
%       names will be 's1205-03-d001.mat','s1205-03-d013.mat' and
%       's1205-03-d114.mat'.
%
%   See also TZ_LOADFEATS_DS

%   ??-???-???? Initial write TINGZ
%   25-MAR-2003 Modified TINGZ
%   14-DEC-2003 Modified TINGZ
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argumentis required')
end

feats=[];
 
if length(label) > 0
   
    labelnum=length(label)
    for j=1:labelnum
        filenum=length(label{j}.num)
        for i = 1:filenum
            labelnum='000'
            no=num2str(label{j}.num(i));
            no
            labelnum(end-length(no)+1:end)=no;
            fullpathname=[dirname '/' label{j}.prefix labelnum '.mat']
            load(fullpathname)
            feats=[feats image_features'];
        end
    end 
else
    feat_list=ml_dir(dirname);
    feat_list=feat_list(3:end);
    feat_num=length(feat_list);
    for i=1:feat_num
        if feat_list{i}(1:2) ~='ts'
          fullpathname=[dirname '/' feat_list{i}]
          load(fullpathname)
          feats=[feats image_features'];
        end
    end
end

feats=feats';
