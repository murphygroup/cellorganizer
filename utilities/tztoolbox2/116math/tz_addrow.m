function B=tz_addrow(A,x)
%TZ_ADDROW Obsolete.

%function B=tz_subrow(A,x)

error('tz_addrow is out of date. Please use ml_addrow');

B=A+ones(size(A,1),1)*x;