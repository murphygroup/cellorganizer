function feat = tz_medaxspfeat(medaxis,width,param)
%TZ_MEDAXSPFEAT Obsolete. See ML_MEDAXSPFEAT.
%   FEAT = TZ_MEDAXSPFEAT(MEDAXIS,WIDTH) returns a cell array that has 5
%   elements.. The first element is the length of the medial axis. FEAT{2}
%   and FEAT{3} contain the internal nodes and coeffiecients for spline fit 
%   with the medial axis MEDAXIS. FEAT{4} and FEAT{5} contain the internal 
%   nodes and coeffiecients for for spline fit with width WIDTH. Here both
%   curves have one internal node and 5 coefficeints.
%   
%   FEAT = TZ_MEDAXSPFEAT(MEDAXIS,WIDTH,PARAM) specifies the spline to fit
%   with MEDAXIS and WIDTH. PARAM is a structure and it has two fields,
%   knot sequence 'knots' and spline order 'k'. See SPAP2 for more details.
%   
%   See also

%   27-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

error(tz_genmsg('of','tz_medaxspfeat','ml_medaxspfeat'));

if nargin < 2
    error('2 or 3 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param, struct('knots1',2,...
    'k1',4,'knots2',2,'k2',4,'isshow',0,'w',ones(1,size(medaxis,1))));

y = medaxis(:,2)';
len = length(y);
feat{1} = len;
x=(0:(len-1))/(len-1);

%flip the data upside down if necessary
y1=y(x<0.33 | x>=0.67);
y2=y(x>=0.33 & x<0.67);
if mean(y1)>mean(y2)
    y=-y;
end
y=y-min(y);

sp=spap2(param.knots1,param.k1,x,y,param.w);
[feat{2},feat{3}] = tz_sp2feat(sp);
    
if param.isshow
    subplot(1,2,1)  
    plot(x,y,'rx');
    hold on
    fnplt(sp,'k-',2);
    hold off
end

y = width;

sp=spap2(param.knots2,param.k2,x,y,param.w);
[feat{4},feat{5}] = tz_sp2feat(sp);
   
if param.isshow
    subplot(1,2,2)
    plot(x,y,'rx');
    hold on
    fnplt(sp,'k-',2);
    hold off
end
