function results = tz_classif2_df(featsfile,task)
%TZ_CLASSIF2_DF Binary classification for feautres of drug images.
%   RESULTS = TZ_CLASSIF2_DF(FEATSFILE,TASK) returns the average accuracy
%   of classifying group pairs from the feature file FEATSFILE. TASK
%   specifys how the group pairs are formed:
%       'between' - control group and drug group
%       '~within' - two parts of control group
%       'within' - two parts of drug group
%       'random' - two parts of mixture of control group and drug group
%
%   RESULTS = TZ_CLASSIF2_DF(FEATSFILE) uses the TASK 'between'.

%   ??-???-???? Initial write T. Zhao
%   27-MAR-2003 Modified T. Zhao
%   24-DEC-2003 Modified T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University


if ~exist('task','var')
     task = 'between';
end

load(featsfile)

nprot = size(features,3);
ndrug = size(features,2);
m=1;
for i=1:nprot
    for j=1:ndrug
        switch(task)
            case 'between',
                f1 = features{1,j,i};
                f2 = features{2,j,i};
            case '~within',
                samplenum=size(features{1,j,i},1);
                partlen=round(samplenum/2);
                ranorder=tz_randorder(samplenum);
                f1=features{1,j,i}(1:partlen,:);
                f2=features{1,j,i}(partlen+1:end,:);
            case 'within',
                samplenum=size(features{2,j,i},1);
                partlen=round(samplenum/2);
                ranorder=tz_randorder(samplenum);
                f1=features{2,j,i}(1:partlen,:);
                f2=features{2,j,i}(partlen+1:end,:);
            case 'random'
                f=[features{1,j,i}; features{2,j,i}];
                samplenum=size(f,1);
                partlen=samplenum/2;
                ranorder=tz_randorder(samplenum);
                f1=f(uint16(ranorder(1:round(partlen))),:);
                f2=f(uint16(ranorder(round(partlen)+1:end)),:);
                %round(partlen)
            otherwise,
                'incorrect task name';
        end
        
        
        if( ~isempty(f1) & ~isempty(f2))
            tz_make_classfile('tmpallfeats.mat',{f1,f2});
            correct=hk_exprbfsvm('tmpallfeats.mat');
            delete tmpallfeats.mat
            sigcor(m)=correct;
            m=m+1;
        else
            correct = nan;
        end
        results(j,i) = correct;
    end
end

mean_correct=mean(sigcor)