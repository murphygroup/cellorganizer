function regmodel = tz_ldareg(x,y,t)
%TZ_LDAREG Obsolete.

%function regomdel = tz_ldareg(x,y,t)
%OVERVIEW
%   regression by lda
%PARAMETERS
%   x - variables
%   y - targets
%   t - regression parameters
%RETURN
%   regmodel - regression model
%DESCRIPTION
%   Actually this is only for classifying categorized data.
%HISTORY
%   29-May-2005 Initial write TINGZ
%SEE ALSO
%   tz_evalldareg

error(tz_genmsg('of','tz_ldareg','ml_ldareg'));

if (size(y,2)>1)
    error('The target y must have one column');
end

ngroup = max(y);
gr=size(y,1);

constidx=[];
for i=1:ngroup
    gx=x(y==i,:);
    constidx=[constidx tz_constidx(gx)];
end

if isempty(constidx)
    prep.featidx=[];
else
    allidx=1:size(x,2);
    allidx(constidx)=[];
    prep.featidx=allidx;
    x(:,constidx)=[];
end
    
[rx,cx] = size(x);
t.trsize=rx;


if rx ~= gr,
    regmodel.succ=0;
    regmodel.warnmsg='The number of rows in the second and third input arguments must match.';
    error(regmodel.warnmsg); 
end

for k = 1:ngroup
   groupk = x(find(y == k),:);
   [rg,cg]=size(groupk);
   if rg < cg
       regmodel.succ=0;
       regmodel.warnmsg='The number of samples must exceed the number features for LDA.';
       error(regmodel.warnmsg); 
   end
   meanx = mean(groupk);
   [Q,R] = qr(groupk - meanx(ones(rg,1),:),0);
   trained.rs{k}=R;
   trained.means{k}=meanx;
   trained.rx(k)=rg;
end

regmodel.modelname='lda';
regmodel.type='llk'; %likelihood
regmodel.trained=trained;
regmodel.t=t;
regmodel.prep=prep;
regmodel.postp.ctg=1;

