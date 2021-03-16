function regmodel = ml_regress(x,y,t,regmeth)
%ML_REGRESS Multivariate regression.
%   REGMODEL=ML_REGRESS(X,Y,T,REGMETH) returns a structure restoring 
%   the regression model of f(X)=E(Y|X). T is the structure of
%   parameters. It has the following fields:
%       norm: 0 for no normalization,
%             1 for normalization upon all training data,
%             2 for normalization upon training set
%               if the training data is splitted into a training
%               set and a stop set
%       stop: 0 for no stop set
%             1 for using stop set. 
%       randtrainsel (optional): indices of training set in training data
%   If T is empty, default parameters will be used.
%   REGMETH is an integer specifying regression method:
%       1 for Linear Discriminant Analysis (LDA) (See ML_LDAREG)
%       2 for Support Vector Machine (SVM) (See ML_SVMREG)
%       3 for Back Propagation Neural Network (BPNN) (See ML_BPNNREG)
%       4 for LIBSVM (See ML_LIBSVMREG)
%
%   SEE ALSO ML_EVALREG

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

if isempty(t)
    t = struct([]);
end
    
switch regmeth
    case 1 %LDA
        regmodel = ml_ldareg(x,y);
        regmodel.t = t;
    case 2 %SVM
        t = ml_initparam(t,struct('norm',1,'stop',0,'C_values',20, ...
            'rbf_levels',7,'model_types',maxwin,'tutor',smosvctutor));
        regmodel=ml_svmreg(x,y,t);
    case 3 %BPNN
        t = ml_initparam(t,struct('norm',2,'stop',1,'epochs',300, ...
            'bp',1000, 'hidden',20));
        regmodel = ml_bpnnreg(x,y,t);
    case 4 %LIBSVM
        t = ml_initparam(t, ...
              struct('norm',3,'args','-s 0 -b 1 -c 512.0 -g 0.0078125'));
        regmodel = ml_libsvmreg(x,y,t);
    otherwise
        error(['Unrecognized regression option: ' num2str(regmeth)]);
end

