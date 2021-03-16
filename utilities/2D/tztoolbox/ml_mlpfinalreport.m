function ml_mlpfinalreport(trainnetout,trainclass,trainidx, ...
			   stopnetout,stopclass,stopidx,testnetout,...
			   testclass,testidx,classnames,filenames, fileroot)
%  ML_MLPFINALREPORT -
%
%  [] = ML_MLPFINALREPORT()
%
%    Outputs:
%
%    Inputs:
%
%
%    M. Boland - 21 May 1999
%

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

% $Id: ml_mlpfinalreport.m,v 1.2 2006/06/27 13:33:47 tingz Exp $

%
% Find the best thresholds for each network trial
%
thresholds = [0.05:0.05:0.95] ;
onlyone=1 ;
[bestthresh] = ml_mlpthreshtest(stopnetout,stopclass,thresholds,onlyone) ;
%
% Combine the confusion matrices from the network trials
%
[trainconfmats]=ml_mlpconfmatall(trainnetout,trainclass,bestthresh,onlyone) ;
[stopconfmats] = ml_mlpconfmatall(stopnetout,stopclass,bestthresh,onlyone) ;
[testconfmats] = ml_mlpconfmatall(testnetout,testclass,bestthresh,onlyone) ;
%
% Summarize the confusion matrices
%
[trainsummary] = ml_mlpclasssummary(trainconfmats) ;
[stopsummary] = ml_mlpclasssummary(stopconfmats) ;
[testsummary] = ml_mlpclasssummary(testconfmats) ;
%
% write the output to a file
%
diary(strcat(fileroot,'_conf_thresh.txt')) ;
disp('Training Data')
trainsummary.confusion 
trainsummary.confusion_nounk
trainsummary.Pc_mean
trainsummary.Pc_var

disp('Stop Data')
stopsummary.confusion 
stopsummary.confusion_nounk
stopsummary.Pc_mean
stopsummary.Pc_var

disp('Test Data')
testsummary.confusion 
testsummary.confusion_nounk
testsummary.Pc_mean
testsummary.Pc_var

diary
%
% write LaTeX tables for the test data
%

 classnames = {'DNA','ER','Giantin','GPP130','LAMP2','Mito.','Nucle.','Actin','TfR','Tubul.'}
ml_latextable(strcat(fileroot,'_test_thresh.tex'),testsummary.confusion*100,...
              classnames) ;
ml_latextable(strcat(fileroot,'_test_thresh_nounk.tex'),...
              testsummary.confusion_nounk*100, classnames) ;

%
% Combine the confusion matrices from the network trials, NO thresholding
%
 
[trainconfmats] = ml_mlpconfmatall(trainnetout,trainclass,0,0) ;
[stopconfmats] = ml_mlpconfmatall(stopnetout,stopclass,0,0) ;
[testconfmats] = ml_mlpconfmatall(testnetout,testclass,0,0) ;
%
% Summarize the confusion matrices
%
[trainsummary] = ml_mlpclasssummary(trainconfmats) ;
[stopsummary] = ml_mlpclasssummary(stopconfmats) ;
[testsummary] = ml_mlpclasssummary(testconfmats) ;
%
% write the output to a file
%
diary(strcat(fileroot,'_conf.txt')) ;
disp('Training Data')
trainsummary.confusion 
trainsummary.confusion_nounk
trainsummary.Pc_mean
trainsummary.Pc_var

disp('Stop Data')
stopsummary.confusion 
stopsummary.confusion_nounk
stopsummary.Pc_mean
stopsummary.Pc_var

disp('Test Data')
testsummary.confusion 
testsummary.confusion_nounk
testsummary.Pc_mean
testsummary.Pc_var

diary
%
% write LaTeX tables for the test data
%
ml_latextable(strcat(fileroot,'_test.tex'),testsummary.confusion*100, ...
	      classnames) ;
ml_latextable(strcat(fileroot,'_test_nounk.tex'),...
              testsummary.confusion_nounk*100, classnames) ;




[testnames missed assignedto] = ml_missclassnames(filenames,testidx,testconfmats) ;

ml_missclassnamesprint(strcat(fileroot,'_missed.txt'),missed,assignedto) ;

