function params = tz_trainbayesmnom(x,isshow)
%TZ_TRAINBAYESMNOM Train a bayesian multinomial distribution.
%   PARAMS = TZ_TRAINBAYESMNOM(X) returns a structure containing parameters
%   for the mutimonmial distribution with Dirichlet distribution from 
%   empirical Bayes estimation.
%
%   params = tz_trainbayesmnom(x,1) will also show real-time training
%   likelihood.
%   
%   See also

%   21-Oct-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

if nargin<2
    isshow = 0;
end

x = sum(x,1);
k = size(x,2);

%initialize
RHO = 0.1;
MINERROR = 1e-5;
MAXITER = 1000;
totalNumber = sum(x);
lamdas = ones(1,k);
lamdas(1)=5;
oldLamdas = lamdas+MINERROR*2;
niter = 1;
density = [];
% cx = factorial(sum(x))/prod(gamma(x+1));
while any(abs(oldLamdas-lamdas)>MINERROR)
    oldLamdas = lamdas;
    oldAlphas = exp(oldLamdas);
    oldA = sum(oldAlphas);
    sumLamda = sum(oldLamdas);
    
    if isshow
        density = [density dcnormfactor(oldAlphas)/ ...
                dcnormfactor(oldAlphas+x)];
        plot(density);  
        drawnow
    end

    lamdas = oldLamdas;
    alphas = oldAlphas;
    for i=1:k
        lamdas(i) = lamdas(i) + RHO*exp(lamdas(i)).*(digamma(oldA)- ...
            digamma(oldA+totalNumber)- ...
            digamma(oldAlphas(i))+ ...
            digamma(oldAlphas(i)+x(i)));
        alphas = exp(lamdas);
        oldA = sum(alphas);
    end

    niter = niter+1;

    if niter>MAXITER
        warning('Maximum iteration exceeded.');
        break;
    end
end
% alphas = ones(1,k);
% while any(abs(oldLamdas-lamdas)>MINERROR)
%     oldAlphas = alphas;
%     oldA = sum(oldAlphas);
%     
%     if isshow
%         density = [density dcnormfactor(x)*dcnormfactor(oldAlphas)/ ...
%                 dcnormfactor(oldAlphas+x)];
%         plot(density);  
%         drawnow
%     end
%     (digamma(oldA)- ...
%         digamma(oldA+totalNumber)- ...
%         digamma(oldAlphas)+ ...
%         digamma(oldAlphas+x))
%     alphas = oldAlphas + RHO.*(digamma(oldA)- ...
%         digamma(oldA+totalNumber)- ...
%         digamma(oldAlphas)+ ...
%         digamma(oldAlphas+x));
%  
%     niter = niter+1;
%     if niter>MAXITER
%         warning('Maximum iteration exceeded.');
%         break;
%     end
% end

alphas = exp(lamdas);
theta = (x+alphas)/(sum(alphas)+totalNumber);
params = struct('alpha',alphas,'theta',theta);

function c=dcnormfactor(alpha)
%function dcnormfactor(alpha)
c = gamma(sum(alpha))/prod(gamma(alpha));
    