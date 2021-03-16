function p_value = tp_pvalue(x,f)
% Calculate the p-value if a random variable(vector) X is drawn from a
% probability distribution F
%Edited: 
%7/6/13 - D.Sullivan added catch for if variance==0

if nargin < 2
    error('At least 2 argument is required');
end

%D. Sullivan 7/6/13 - added support for sigma==0
if all(f.sigma(:) == 0) && all(x ==f.mu)
    p_value = 1;
    return
elseif all(f.sigma(:) == 0)
    p_value = 0;
    return
end

switch f.name
    case 'norm'
        xcdf = normcdf(x,f.mu,f.sigma);
        p_value = 2*xcdf*(xcdf<.5)+2*(1-xcdf)*(xcdf>=.5);
    case 'mvn'
        x = x - f.mu;
        chi_stat = x * (f.sigma\x');
        p_value = 1 - chi2cdf(chi_stat,length(x));
%     case 'gamma'
%     case 'exp'
%     case 'mix'
    otherwise
        error('Unrecognized pdf.');
end
