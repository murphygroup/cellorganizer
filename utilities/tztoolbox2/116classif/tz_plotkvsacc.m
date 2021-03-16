function tz_plotkvsacc(datadir,set,k)
%TZ_PLOTKVSACC Plot accuracies vs numbers of object types.
%   TZ_PLOTKVSACC(DATADIR,SET,K) shows the graph of classification
%   accuracy vs number of object types. The plot data are read from
%   DATADIR. SET is a string of classifiction infomation, which is used
%   to find proper data files. K is a vector of number of object types.

%   02-APR-2004 Initial write T. Zhao
%   05-NOV-2004 Modified T. Zhao
%       - add comments
%       - change parameter dir to datadir
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 3
    error('Exactly 3 arguments are required')
end

for i=1:length(k)
    load([datadir '/' 'results' num2str(k(i)) set '10fold.mat']);
    [avgcm2{i},tempavg]=tz_calcavgcm(cvcm2);
    avgacc(i)=mean(avgacc2);
        
end
[maxacc,pos]=max(avgacc);
plot(k,avgacc,'.');
xlabel('k')
ylabel('average accuracy')
hold on
plot(k(pos),maxacc,'r+');
hold off