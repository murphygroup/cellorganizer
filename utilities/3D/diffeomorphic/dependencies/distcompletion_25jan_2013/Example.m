% Test script to run low-rank matrix completion algorithm from
% B. Mishra, G. Meyer and R. Sepulchre, "Low-rank optimization for distance matrix completion", 
% Submitted to the 50th IEEE Conference on Decision and Control, December 2011

% Last modified, april 2011

clear all; close all; clc;

% Set randstream based on clock
s = RandStream.create('mt19937ar','seed',sum(100*clock));
RandStream.setGlobalStream(s);

% Select optimization algorithm
%methodName = 'Gradient Descent';
methodName = 'Trust Region';   

n = 500;

% Fraction of unknown distances
fractionOfUnknown = 0.80;

% Rank of the solution
r = 5;
% Starting approximation rank
p = 1;

% Parameters
params.pmax = r;
params.tol = 1e-3;
params.vtol = 1e-3;
params.verb = true; % "Bark" only when asked

switch methodName,
  case 'Gradient Descent',
    methodFun = @gd_dist_completion;
  case 'Trust Region',
    methodFun = @tr_dist_completion;
  otherwise,
    error('Unknown method %s',methodName);
end

Yo = randn(n,r); % true embedding
trueDists = pdist(Yo)'.^2; % true distances

% Comment out this line if you don't want to add noise
trueDists = trueDists + 0.01 * std(trueDists) * randn(size(trueDists)); % add noise

% Compute all pair of indices
H = tril(true(n),-1);
[I,J] = ind2sub([n,n],find(H(:))); 
clear 'H';

% HoldOut train-test split
test = false(length(trueDists),1);
test(1:floor(length(trueDists)*fractionOfUnknown)) = true;
test = test(randperm(length(test)));
train = ~test;

% Initial condition
Y0 = randn(n,p);

% Run algorithms
t0 = tic;
[Y infos] = lowrank_dist_completion(methodFun,I(train),J(train),trueDists(train),Y0,params);
timeTaken = toc(t0);

nIter = length(infos.costs) - 1;
        
% Compute predictions
estimDists = pdist(Y)'.^2;

% Compute mean squared errors (MSE)
errorOnTest = mean((trueDists(test) - estimDists(test)).^2);
errorOnTrain = mean((trueDists(train) - estimDists(train)).^2);
errorOverall = mean((trueDists - estimDists).^2);

fprintf('Time taken = %f sec.\n',timeTaken);
fprintf('MSE on test set = %f\n',errorOnTest);
fprintf('MSE on train set = %f\n',errorOnTrain);
fprintf('Total MSE (train + test) = %f\n',errorOverall);

figure;
plot(infos.costs);
legend(sprintf('%s',methodName));
xlabel('Number of iterations');
ylabel('Cost function');
