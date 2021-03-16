function [rnetout,tnetout,imin,net] = tz_mlptraintest(ri,ro,ti,to,hidden,epochs,dp)
%TZ_MLPTRAINTEST Obsolete.

%function [rnetout,tnetout,imin,net] = tz_mlptraintest(ri,ro,ti,to,hidden,epochs,dp)
%OVERVIEW
%   
%PARAMETERS
%   ri - 
%   ro - 
%   ti - 
%   to - 
%   hidden - 
%   epochs - 
%   dp - 
%RETURN
%   rnetout - 
%   tnetout - 
%   imin - 
%   net - 
%DESCRIPTION
%   It does the same job as ml_mlptraintest, but it is more flexible by taking more papameters.
%HISTORY
%   15-May-2005 Initial write TINGZ
%SEE ALSO
%   ml_mlptraintest
error(tz_genmsg('of','tz_mlptraintest','ml_mlptrainteset'));

net = mlp(size(ri,2), hidden, size(ro,2), 'logistic') ;

roptions = zeros(1,18) ;
roptions(1) = 1 ;   % Output sse values
%roptions(1) = -1 ;  % Output nothing 
roptions(14) = 1 ;  % Number of epochs (train one epoch at a time)
roptions(17) = 0.9 ;
roptions(18) = 0.001 ;

[net,rsse,tsse] = ml_mlptrain(net, roptions, ri, ro, ti, to, epochs) ;

%
% Round tsse to the nearest 0.001 to avoid long training sessions
%  with little progress
%
if ~exist('dp','var')
    dp=1000;
end

[min,imin] = min(round(tsse*dp)/dp) ;
pass=1 ;

%
% Continue training until the minimum SSE occurs less than 90% 
%  of the way through the last pass.
%
while ((imin ./ pass) > 0.9*epochs)
  rssesave = rsse ;
  tssesave = tsse ;
  [net,rsse,tsse] = ml_mlptrain(net,roptions,ri,ro,ti,to,epochs) ;
  rsse = [rssesave rsse] ;
  tsse = [tssesave tsse] ;
  [min,imin] = min(round(tsse*dp)/dp) ;
  drawnow
%   plot([1:length(rsse)], rsse, 'r-', [1:length(tsse)], tsse, 'b-');
  pass=pass+1 ;
  pass
  min
  imin
end

%plot([1:length(rsse)], rsse, 'r-', [1:length(tsse)], tsse, 'b-') ;

net = mlp(size(ri,2), hidden, size(ro,2), 'logistic') ;

[net,rsse,tsse] = ml_mlptrain(net, roptions, ri, ro, ti, to, imin) ;

rnetout = mlpfwd(net, ri) ;
tnetout = mlpfwd(net, ti) ;

