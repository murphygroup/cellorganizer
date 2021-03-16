function sig=tz_evalmixsig(ealpha,talpha,B,sigcalc)
%TZ_EVALMIXSIG Evaluate the significance of mixture model estimation.
%   SIG = TZ_EVALMIXSIG(EALPHA,TALPHA,B) returns the significance level of
%   the mixture decomposition results. EALPHA is the estimated component
%   weights and TALPHA is the ground truth. B is the number of trials to
%   estimate the significance level.
%   
%   SIG = TZ_EVALMIXSIG(EALPHA,TALPHA,B,SIGCALC) specifies how to calculate
%   the significance level by SIGCALC. If SIGCALC is 'norm', which is default,
%   the significance level is calculated based on a normal distribution. If it
%   is nonp, the value is calculated based on percentages.
%   
%   See also

%   26-MAY-2004 Initial write TINGZ
%   Copyright (c) Murphy Lab, Carnegie Mellon University


if ~exist('sigcalc','var')
    sigcalc='norm';
end

rbias=tz_simmixbias(talpha,B);

bias=sum(abs(ealpha-talpha));

switch sigcalc
case 'norm'
    sig=normcdf(bias,mean(rbias),sqrt(var(rbias)));
case 'nonp'
    sig=sum(rbias<bias)/B;
end

