function [ratio,interval]=tz_oddsratio(data,sig)

%function [ratio,interval]=tz_oddsratio(data,sig)
%
%OVERVIEW:
%   calculate odds ratio of 2x2 table
%PRAMETERS:
%   data - 2x2 table
%   sig - significance of confidence interval
%RETURN:
%   ratio - odds ratio
%   interval - confidence interval

ratio=(data(1,1)*data(2,2))/(data(1,2)*data(2,1));
lratio=log(ratio);
bias=norminv(1-sig/2,0,1);
sigma=sqrt(sum(1./data(:)));
linterval=[lratio-bias*sigma,lratio+bias*sigma];
interval=exp(linterval);