function [pvalues, tss, boths, BHboths, BHpvalues]= ...
    tz_compfeats_df(features,names,siglevel,task,utestmeth,mtestmeth,t)
%TZ_COMPFEATS_DF Compare all features in df.
%   PVALUES = TZ_COMPFEATS_DF(FEATURES,NAMES,SIGLEVEL) return p-values
%   of pairs of groups in FEATURES. NAMES are protein and drug names.
%   SIGLEVEL is the significance level. The multivariate test is the
%   Hotelling T2 test and the univariate test is the t test.
%   
%   PVALUES = TZ_COMPFEATS_DF(FEATURES,NAMES,SIGLEVEL,TASK) specifies
%   how to form a pair of group from FEATURES by TASK:
%       'between' - control group and drug group
%       '~within' - random groups splitted from control group
%       'within' - random groups splitted from drug group
%       'random' - random groups from control group and drug group
%                   for the same protein
%   
%   PVALUES = TZ_COMPFEATS_DF(FEATURES,NAMES,SIGLEVEL,TASK,UTESTMETH)
%   lets users specify the univariate testing method. See 
%   ML_COMPFEATS_U for more details.
%   
%   PVALUES = TZ_COMPFEATS_DF(FEATURES,NAMES,SIGLEVEL,TASK,UTESTMETH,
%   MTESTMETH,T) lets users specify the multivariate testing method
%   and corresponding parameters. See ML_COMPFEATS_M for more details.
%   
%   [PVALUES,TSS,BOTHS,BHBOTHS,BHPVALUES] = TZ_COMPFEATS_DF(...)
%   also returns test statistics, univariate testing results, 
%   results of BH method. See ML_MULTEST2 for more details.

%   ??-???-???? Initial write T. Zhao
%   27-MAR-2003 Modified T. Zhao
%   09-JUL-2003 Modified T. Zhao
%   14-JUL-2003 Modified T. Zhao
%   14-dec-2003 Modified T. Zhao
%   13-JAN-2003 Modified T. Zhao
%   21-SEP-2003 Modified T. Zhao
%       - add multivariate test methods
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('At least 3 arguments are required')
end

if ~exist('task','var')
     task = 'between';
end

if ~exist('utestmeth','var')
    utestmeth='ml_ttest2';
end

if ~exist('mtestmeth','var')
    mtestmeth='ml_ht2test2';
    t={};
end

nprot = size(features,3);
ndrug = size(features,2);

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
        
        if(length(names)>0)
            names{2,j,i}
        end
        
        if( ~isempty(f1) & ~isempty(f2))
            [pvalue,ts]=ml_compfeats_m(f1,f2,mtestmeth,t)
            [size(f1,1), size(f2,1)]
            %[F, crit, res, both, succ, hotel_warn] = er_compfeat(f1,f2,siglevel);
            %pvalue=tz_ht2test2(f1,f2,0);
            both=ml_compfeats_u(f1,f2,utestmeth,0);
            [BHpvalue,BHboth]=ml_multest2(f1,f2,'BH','ready',0,siglevel,both(2,:));
            %BHresult=any(BHboth(2,:)==1);
        else
            F = nan;
            pvalue=nan;
            both=nan;
            BHboth=nan;
            BHpvalue=nan;
        end 
        
        %HTresults(j,i) = res;
        tss(j,i) = ts;
        pvalues(j,i)=pvalue;
        boths{j,i}=both;
        BHboths{j,i}=BHboth;
        BHpvalues(j,i)=BHpvalue;
    end
end

% num_ones = length(find(HTresults==1));
