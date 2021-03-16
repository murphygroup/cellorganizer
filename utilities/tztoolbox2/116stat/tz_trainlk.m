function params=tz_trainlk(X,trainparam)
%TZ_TRAINLK Train a statistical model for calculating likelihood for data
%   PARAMS = TZ_TRAINLK(X,TRAINPARAM) returns a structure that contains
%   parameters trained from data X. TRAINPARAM is a structure which
%   specifies the parameters for training. It has a field 'method' which
%   spcifies the statistical model:
%       'mnom' - multinomial      
%       'mnomn' - multinomial with numbers
%       'mvn' - normal distribution
%       'mvnw' - weighted normal distribution
%           'weights' - 
%           't' - 
%       'mvnm' - gaussian mixture
%       'mvno' - object type gaussian mixture. The last column of X will be
%           object type [label vector]. 
%       'pca' - principal component analysis
%       'exp' - exponential distribution (only for univariate)
%
%   The field 't' specifies additional parameters , which is a cell array. 
%   Currently it is only available for 'mvnw' and 'mvno' to be the weights.
%  
%   See also TZ_EVALK

%   20-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if ~exist('trainparam','var')
    trainparam = struct([]);
end

trainparam = ml_initparam(trainparam, ...
    struct('method','mvn','weights',[],'t',[]));

% if ~isfield(trainparam,'weights')
%     trainparam.weights = [];
% end

if length(trainparam.method)>=4
    if strcmp(trainparam.method(1:4),'mnom')
        totalnum=sum(X,1);
        p=totalnum/sum(totalnum);
    end
    
    if strcmp(trainparam.method(1:4),'mvno')
        objtypes=X(:,end);
        X(:,end)=[];   
    end
    
    if strmatch(trainparam.method(1:4),{'mvno','mvnw'}) 
        t = trainparam.t;
        if ~isempty(trainparam.weights)            
            weights = t{1}{t{2}};
        else
            weights = [];
        end
    end
end

if strcmp(trainparam.method(1:3),'mvn')
    constidx=ml_constidx(X);
    constx=X(1,constidx);
    tX=X;

    if ~isempty(constidx)
        tX(:,constidx)=[];
    end
end

switch(trainparam.method)
    case 'mnom'
        params=struct('method',trainparam.method,'p',p);
    case 'mnomn'    %mutinomial with numbers
        objnum=[sum(X,2),ones(size(X,1),1)];
        params=struct('method',trainparam.method,'p',p,'objnum',objnum);
    case 'mvn'      %multivariate normal distribution
        mu=mean(tX,1);
        sigma=cov(tX);

        params=struct('method',trainparam.method,'constidx',constidx,'mu',mu, ...
            'sigma',sigma,'constx',constx);

    case 'mvnw'
        if ~isempty(t)
            if ~iscell(t{1})
                weights=t{1};
            else
                weights=t{1}{t{2}};
            end
        else
            weights = ones(size(tX,1),1);
        end

        mu=tz_mean(tX,weights);
        sigma=tz_cov(tX,0,weights);
        params=struct('method',trainparam.method,'constidx', constidx, ...
            'mu',mu, 'sigma',sigma,'constx',constx);

    case 'mvnm' %gaussian mixture
        [bestk,bestpp,bestmu,bestcov,dl,countf] = mixtures4(tX',1,4,0.01,1e-6,0,0);

        params=struct('method',trainparam.method,'constidx', ...
            constidx,'m',bestk,'mus',bestmu,...
            'covs',bestcov,'pp',bestpp,'constx',constx);

    case 'mvno' %obj type gaussian mixture
        trainparam = ml_initparam(trainparam,struct('tz_estcov', ...
            struct('method','mle')));
        
        m=max(objtypes);

        for i=1:m
            pp(i)=sum(objtypes==i);
            if pp(i)~=0
                sx=tX(objtypes==i,:);
                if ~isempty(weights)
                    sweights = weights(objtypes==i,:);
                    pp(i) = pp(i)*mean(sweights);
                end
                mus(:,i)=mean(sx,1)';

                if(size(sx,1)==1) %give perturbation to groups with 1 sample
                    covs(:,:,i)=eye(length(sx))*abs(min(sx(sx~=0)))/100;
                else
                    covs(:,:,i)=tz_estcov(sx,trainparam.tz_estcov);
                end
            end

        end
        mus(:,pp==0)=[];
        covs(:,:,pp==0)=[];
        pp(pp==0)=[];

        pp=pp/sum(pp);
        params=struct('method',trainparam.method,'constidx', ...
            constidx,'m',length(pp),'mus',mus,...
            'covs',covs,'pp',pp,'constx',constx);
    
    case 'pca'  %principal component analysis
        trainparam = ml_initparam(trainparam,struct('ncomp',size(X,2)));
        params.method = trainparam.method;
        params.mean = mean(X,1);
        [pcavec,coeff] = princomp(X);
        params.pcavec = pcavec(:,1:trainparam.ncomp)';
        coeff = coeff(:,1:trainparam.ncomp);
        params.coeffcov = cov(coeff);
        params.coeffmean = mean(coeff,1);
    case 'exp' %exponential distribution
%         beta = tz_fitexp(X,0);
        beta = expfit(X);
        params = struct('method',trainparam.method,'beta',beta);
    case 'poiss' %poisson distribution
        lamda = poissfit(X);
        params = struct('method',trainparam.method,'lamda',lamda);
    case 'norm' %univariate normal distribution
        [params.mu,params.sigma] = normfit(X);
        params.method = trainparam.method;
    otherwise
        error('Unrecognized lk method');
end
