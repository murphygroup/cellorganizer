function regmodel = ml_trainclassif(x,y,t,regmeth)
%ML_TRAINCLASSIF train a classifier
%   REGMODEL=ML_TRAINCLASSIF(X,Y,T,REGMETH) trains a classifier. The returned
%   classifier REGMODEL can be used for the function ML_EVALREG.
%   The inputs are very like those in for the function ML_REGRESS. But T has
%   a field 'featsel' for feature selection. If the value of 'featsel' is 
%   'none', no feature selection will be done. If it is 'sda', the stepwise
%   discriminant analysis (SDA) will be used to select features in each fold.
%   If could also be a vector, which contains the indices of features to
%   select.
%
%   SEE ALSO ML_EVALREG ML_REGRESS

% Copyright (C) 2006  Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

%   12-Jul-2005 Initial write TINGZ
%   26-SEP-2005 Modified TINGZ
%       - just add a semicolon

if max(y) >= 2
    if sum(y==1) == sum(y==2)
        if all( x(y==1,:) == x(y==2,:) )
            nonzeroIndices = find(x~=0);
            if isempty(nonzeroIndices)
                x(1) = x(1)+0.01;
            else
                x( nonzeroIndices(1) ) = x( nonzeroIndices(1) ) * 1.01;
            end
        end
    end
end

t = ml_initparam(t,struct('featsel','none'));
prep.featidx = [];
switch num2str(t.featsel)
    case 'sda'
        prep.featidx = ml_stepdisc(ml_combfeats2mcf(x,y));
    otherwise
        if isnumeric(t.featsel)
            prep.featidx = t.featsel;
        end
end

allidx = 1:size(x,2);
if ~isempty(prep.featidx)
    x = x(:,prep.featidx);
    allidx = allidx(prep.featidx);
end


switch regmeth
    case 1
%     regmodel = ml_ldareg(x,y);
    case 2
        y = ml_label2post(y); 
        y(find(y==0)) = -1;  
    case 3
        y=ml_label2post(y)*0.8+0.1;
    case 4
        %do nothing
    otherwise
        error(['Unrecognized  classificatin option: ' num2str(regmeth)]);
end

regmodel = ml_regress(x,y,t,regmeth);
if ~isempty(regmodel.prep.featidx)
    %t-
    %     regmodel.prep.featidx = prep.featidx(regmodel.prep.featidx);
    %t--
    %t+
    regmodel.prep.featidx = allidx(regmodel.prep.featidx);
    %t++
else
    regmodel.prep.featidx = prep.featidx;
end

regmodel.postp.ctg=1;
