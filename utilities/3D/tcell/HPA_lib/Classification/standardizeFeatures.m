function [traindata, testdata, mu, st] = standardizeFeatures( traindata, testdata)
% [TRAINDATA, TESTDATA] = FEATNORM( TRAINDATA, TESTDATA)
% 
% Z-score standardizes data by subtracting column elements by the mean of the 
% respecitve columns, and then dividing by the standard deviation of each column
% 
% Last modified Justin Newberg 11 February 2009

mu = mean(traindata,1);  st = std(traindata,[],1);
traindata = (traindata - repmat( mu, [size(traindata,1) 1]))./ repmat( st, [size(traindata,1) 1]);

if exist('testdata','var')
   testdata = (testdata - repmat( mu, [size(testdata,1) 1]))./ repmat( st, [size(testdata,1) 1]);
else
   testdata  = [];
end
