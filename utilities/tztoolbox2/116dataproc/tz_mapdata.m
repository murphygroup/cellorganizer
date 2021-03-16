function y = tz_mapdata(x,option)
%TZ_MAPDATA Data transformation for texture growing model.
%   Y = TZ_MAPDATA(X,OPTION) returns a vector of transformation of X.
%   This is especially for texture growing model. OPTION is used to
%   specify transformation method. If it is empty, X will be returned
%   without any change.
%   X is usually like:
%   [offset, sequence index, celldist, nucdist, neighbor weights

%   29-Apr-2005  Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

if isempty(option)
    y=x;
    return
end

switch option
case 'gt1'  %[offset, celldist, nucdist]
    y=x;
    y(:,[2,5])=[];
case 'gt2'  %[offset, celldist, nucdist, neighbor weights]
    y=x;
    y(:,2)=[];
case 'gt3'  %[offset, sequence index, celldist, nucdist]
    y=x;
    y(:,5)=[];
case 'gt4'  
    %[offset, sequence index, celldist, nucdist -- isolated
    % offset, sequence index, neighbor weights -- not isolated]
    s1=x(:,5)==0;
    s2=x(:,5)>0;
    y=[x(:,1).*s1,x(:,2).*s1,x(:,3).*s1,x(:,4).*s1,...
            x(:,1).*s2,x(:,2).*s2,x(:,5).*s2];
case 'gt5'
    %[offset, sequence index, celldist, nucdist -- isolated
    % offset, sequence index, celldist, nucdist, neighbor weights -- not isolated]
    s1=x(:,5)==0;
    s2=x(:,5)>0;
    y=[x(:,1).*s1,x(:,2).*s1,x(:,3).*s1,x(:,4).*s1,...
            x(:,1).*s2,x(:,2).*s2,x(:,3).*s2,x(:,4).*s2,x(:,5).*s2];
case 'gt6'
    %[offset, sequence index, celldist, nucdist -- isolated
    % offset, sequence index, neighbor weights, sequence index*neighbor weights-- not isolated]
    s1=x(:,5)==0;
    s2=x(:,5)>0;
    y=[x(:,1).*s1,x(:,2).*s1,x(:,3).*s1,x(:,4).*s1,...
            x(:,1).*s2,x(:,2).*s2,x(:,5).*s2,x(:,5).*x(:,2).*s2];
case 'gt7'
    %[offset, sequence index, celldist, nucdist -- isolated
    % offset, sequence index, neighbor weights, neighbor weights^2-- not isolated]
    s1=x(:,5)==0;
    s2=x(:,5)>0;
    y=[x(:,1).*s1,x(:,2).*s1,x(:,3).*s1,x(:,4).*s1,...
            x(:,1).*s2,x(:,2).*s2,x(:,5).*s2,x(:,5).^2.*s2];
case 'gt8'
    %[offset, sequence index, celldist, nucdist -- isolated
    % celldist*nucdist, celldist ^2, nucdist^2-- isolated
    % offset, sequence index, neighbor weights, neighbor weights^2-- not isolated]
    s1=x(:,5)==0;
    s2=x(:,5)>0;
    y=[x(:,1).*s1,x(:,2).*s1,x(:,3).*s1,x(:,4).*s1,...
            x(:,3).*x(:,4).*s1,x(:,3).^2.*s1,x(:,4).^2.*s1,...
            x(:,1).*s2,x(:,2).*s2,x(:,5).*s2,x(:,5).^2.*s2];
       
case 'gt9'
    %[offset, sequence index, celldist, nucdist -- isolated
    % offset, sequence index, neighbor weights -- not isolated]
    s1=x(:,5)>10;
    s2=x(:,5)<=10;
    y=[x(:,1).*s1,x(:,2).*s1,x(:,3).*s1,x(:,4).*s1,...
            x(:,1).*s2,x(:,2).*s2,x(:,5).*s2];
case 'gt10'
    %[offset, sequence index, celldist, nucdist -- isolated
    % offset, sequence index, celldist, nucdist, neighbor weights -- not isolated]
    s1=x(:,5)>10;
    s2=x(:,5)<=10;
    y=[x(:,1).*s1,x(:,2).*s1,x(:,3).*s1,x(:,4).*s1,...
            x(:,1).*s2,x(:,2).*s2,x(:,3).*s2,x(:,4).*s2,x(:,5).*s2];
case {'gt11','gt11b'}
    %[offset, sequence index, celldist, nucdist -- isolated
    % celldist*nucdist, celldist ^2, nucdist^2 -- isolated
    % offset, sequence index, neighbor weights -- not isolated]
    s1=x(:,5)>10;
    s2=x(:,5)<=10;
    y=[x(:,1).*s1,x(:,2).*s1,x(:,3).*s1,x(:,4).*s1,...
            x(:,3).*x(:,4).*s1,x(:,3).^2.*s1,x(:,4).^2.*s1,...   
            x(:,1).*s2,x(:,2).*s2,x(:,5).*s2];
case 'gt12'
    %[offset, sequence index, celldist, nucdist -- isolated
    % celldist*nucdist, celldist ^2, nucdist^2 -- isolated
    % offset, sequence index, neighbor weights -- not isolated]
    % same as 'gt11'?
    s1=x(:,5)==0;
    s2=x(:,5)>0;
    y=[x(:,1).*s1,x(:,2).*s1,x(:,3).*s1,x(:,4).*s1,...
            x(:,3).*x(:,4).*s1,x(:,3).^2.*s1,x(:,4).^2.*s1,...
            x(:,1).*s2,x(:,2).*s2,x(:,5).*s2,];
case 'gt13'
    %[offset, sequence index, -- isolated
    % offset, sequence index, -- not isolated
    % celldist ]
    y=[x(:,1).*(x(:,end)==0),x(:,2).*(x(:,end)==0),...
            x(:,1).*(x(:,end)>0),x(:,2).*(x(:,end)>0),x(:,3:end)];
 case 'gt14'
     %same as 'gt13'?
    y=[x(:,1).*(x(:,end)==0),x(:,2).*(x(:,end)==0),...
            x(:,1).*(x(:,end)>0),x(:,2).*(x(:,end)>0),x(:,3:end)];
case 'gt15'
    %[offset, sequence index, celldist, nucdist -- isolated
    % offset, sequence index -- not isolated
    % neighbor weights ]
    y=[x(:,1).*(x(:,end)==0),x(:,2).*(x(:,end)==0),x(:,3).*(x(:,end)==0),...
            x(:,4).*(x(:,end)==0),x(:,1).*(x(:,end)>0),x(:,2).*(x(:,end)>0),x(:,end)];
end
