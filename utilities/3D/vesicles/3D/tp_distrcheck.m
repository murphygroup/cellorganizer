function quantile = tp_distrcheck(x,distr)
% Validate if the samples follows a parametric distribution

N = size(x,1);
switch distr
    case 'mvn'
        mvnparam = ml_estpdf(x,struct('name','mvn'));
        for i = 1:size(x,2)
            xaxis = sort(unique(x(:,i)));
            H = hist(x(:,i),xaxis);
            epcq{i} = cumsum(H)/N;
            mu = mvnparam.mu(i);
            sigma = sqrt(mvnparam.sigma(i,i));
            distrq{i} = normcdf(xaxis,mu,sigma);
%             plot(distrq{i},epcq{i},'k-')
%             hold on
%             axis equal, axis([0 1 0 1])
        end
    case 'norm'
        normparam = ml_estpdf(x,struct('name','norm'));
        xaxis = sort(unique(x));
        H = hist(x,xaxis);
        epcq = cumsum(H)/N;
        distrq = normcdf(xaxis,normparam.mu,normparam.sigma);
end

quantile.name = distr;
quantile.epcq = epcq;
quantile.distrq = distrq;