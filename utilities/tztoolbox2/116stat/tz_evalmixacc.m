function acc = tz_evalmixacc(ealpha,talpha)
%TZ_EVALMIXACC Calculate accuracy for mixture decomposition.
%   ACC = TZ_EVALMIXACC(EALPHA,TALPHA)
%   
%   See also

%   19-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function acc = tz_evalmixacc(ealpha,talpha)
%OVERVIEW
%   calculate accuracy for mixture decomposition
%PARAMETERS
%   ealpha - evaluated coefficients
%   talpha - true coefficients
%RETURN
%   acc - accuracy
%DESCRIPTION
%   
%HISTORY
%   16-May-2005 Initial write TINGZ
%SEE ALSO
%   tz_evalmixsig


if any(ealpha<0) | any(ealpha>1) | any(abs(sum(ealpha,2)-1)>1e-5)
    error('invalid proportions');
end
if any(talpha<0) | any(talpha>1) | any(abs(sum(talpha,2)-1)>1e-5)
    error('invalid proportions');
end
a(:,:,1)=ealpha;
a(:,:,2)=talpha;
b=min(a,[],3);
acc=mean(sum(b,2));
