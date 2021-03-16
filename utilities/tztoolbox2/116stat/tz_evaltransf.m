function Y=tz_evaltransf(X,B)

%TZ_EVALTRANSF: evaluate the regression

if size(B,1)-size(X,2)==1
    X=[ones(size(X,1),1),X];
end

Y=X*B;