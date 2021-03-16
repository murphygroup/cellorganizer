function [Y2, R] = mdsClassic(D, dims)
% This function do classical MDS on the input distance matrix D.
% dims is the parameter indicating the dimensions of the reduced dimensionalities

%    Y = Y.coords is a cell array, with coordinates for d-dimensional embeddings
%         in Y.coords{d}.  Y.index contains the indices of the points embedded.
%    R = residual variances for embeddings in Y


D = double(D);
N = size(D,1);
opt.disp = 0; 
verbose = 1;
R = [];
[vec, val] = eigs(-.5*(D.^2 - sum(D.^2)'*ones(1,N)/N - ones(N,1)*sum(D.^2)/N + sum(sum(D.^2))/(N^2)), max(dims), 'LR', opt); 
Y.index = [1:N];

h = real(diag(val)); 
[foo,sorth] = sort(h);  sorth = sorth(end:-1:1); 
val = real(diag(val(sorth,sorth))); 
vec = vec(:,sorth); 

D = reshape(D,N^2,1); 
for di = 1:length(dims)
     if (dims(di)<=N)
         Y.coords{di} = real(vec(:,1:dims(di)).*(ones(N,1)*sqrt(val(1:dims(di)))'))'; 
         Y2.coords{di} = real(vec(:,1:dims(di)))'; 
         r2 = 1-corrcoef(reshape(real(L2_distance(Y.coords{di}, Y.coords{di})),N^2,1),D).^2; 
         R(di) = r2(2,1); 
         if (verbose == 1)
             disp(['  MDS on ' num2str(N) ' points with dimensionality ' num2str(dims(di)) '  --> residual variance = ' num2str(R(di))]); 
         end
     end
end