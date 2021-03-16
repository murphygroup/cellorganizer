function [ d, ss, phis, thetas ] = model2posMap( betas )
%Take the protein location parameter from a cell organizer model
%   1x4 vector corresponding to Beta in T. Zhao and R. F. Murphy (2007). Automated learning of generative models for subcellular location: Building blocks for systems biology. Cytometry 71A:978-990.     
%       OR
%   1x6 vector corresponding to Beta in T. Peng and R.F. Murphy (2011) Image-derived, Three-dimensional Generative Models of Cellular Organization. Cytometry Part A 79A:383-391. 
%
%Returns a 2D or 3D maxtrix of values corresponding to P(s, phi, theta), the
%probability location of an object

%Gregory Johnson late 2012/ early 2013
%G. Johnson 7/7/2013 added support for 2D models

betas = double(betas);

ss = double(-0.5:0.025:1); %r in Zhao 2007
phis = double(-pi:0.1:pi); %a in Peng 2007
thetas = [];

if length(betas) == 6

    thetas = double((-pi/2):0.05:(pi/2));

    [s,t,p] = meshgrid(ss, thetas, phis);
    mappos = ml_mappos([s(:), t(:), p(:)]);

    try
       e = mappos*betas;
    catch
       betas = betas';
       e = mappos*betas;
    end

    %P = probability density associated at each position
    P = exp(e)./(1+exp(e));


else
    P = 0;
    
    rminds = any(isnan(betas),1);
    betas(:,rminds) = [];
    
    [s, p] = meshgrid(ss, phis);
    
    for i = 1:size(betas,2)
        i;
        beta = betas(:,i);
        if ~any(isnan(beta))
            x = [ones(numel(s),1) s(:) s(:).^2 sin(p(:)) cos(p(:))];
            
            pst = ml_evallogistic(x,beta);

            P = P + pst;
        end
    end
    P = P ./ size(betas,2);
end

    d = zeros(size(s));
    d(1:end) = P;


end