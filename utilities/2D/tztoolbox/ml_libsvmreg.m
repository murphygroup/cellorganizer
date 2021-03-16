function regmodel = ml_libsvmreg(x,y,param)
%ML_LIBSVMREG Lib SVM regression.
%   REGMODEL = ML_LIBSVMREG(X,Y) train a classifier by LIBSVM.
%   
%   REGMODEL = ML_LIBSVMREG(X,Y,PARAM) specifies how to train the
%   classifier by PARAM, which is a structure with the following fields:
%       'norm' - data normalization. The default value is 3, which means the
%           training data will be normalized into [-1 1]. Keeping the default
%           value is recommended. 
%       'args' - parameters of the SVM. It is a string containg the arguments
%           of the function SVMTRAIN. The default value is '-s 0 -c 512.0 -g 
%           0.0078125'. To search the optimal parameters automatically, set it
%           to empty ([]).
%
%   Example:
%       x = [rand(50,2); rand(50,2)+0.5];
%       y = [ones(50,1); ones(50,1)+1];
%       regmodel = ml_libsvmreg(x,y,struct('args',[]));
%   
%   See also

%   15-Jan-2007 Initial write T. Zhao
%   Copyright (c) 2007 Murphy Lab
%   Carnegie Mellon University
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation; either version 2 of the License,
%   or (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
%   02110-1301, USA.
%   
%   For additional information visit http://murphylab.web.cmu.edu or
%   send email to murphy@cmu.edu


if nargin < 2
    error('2 or 3 arguments are required');
end

if ~exist('param','var')
    param = struct([]);
end

param = ml_initparam(param, ...
     struct('norm',3,'args','-s 0 -c 512.0 -g 0.0078125'));

[x,prep] = ml_dataprep(x,param);

if isempty(param.args) %search parameters automatically
    prevdir = pwd;
    workdir = [fileparts(which('svmtrain')) '/python'];,
    filepath = tempname;
    feats = ml_combfeats2mcf(x,y);
    ml_makelibsvmfile(filepath,feats);
    cd(workdir);
    [stat,msg]=unix(['python grid.py ' filepath ' | tail -1']);
    %delete *.out
    %delete *.png
    cd(prevdir);
    param.args = parselibsvminfo(msg);  
end

trained = svmtrain(y,x,param.args);

regmodel=struct('modelname','libsvm','modeltype',...
    'svm','trained',trained,'t',param,'prep',prep);

if size(y,2)==1
    regmodel.postp.ctg=1;
end


%%%%%%get svm parameter searching information%%%%%%
function arg = parselibsvminfo(msg)

pos=find(msg==char(10))

if length(pos)>1
    msg = msg(pos(end-1)+1:pos(end)-1);
else
    msg = msg(1:1:pos(end)-1);
end

[bestc,msg] = strtok(msg);
[bestg,msg] = strtok(msg);

arg = ['-s 0 -c ' bestc ' -g ' bestg];
