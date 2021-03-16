function splits = partition( labels, N, rn)
% SPLITS = PARTITION( LABELS, N)
%   a fold partitioning function
%   outputs splits is a cell array with N entries

if ~exist( 'rn','var')
    rn = 13;
end
rand( 'seed',rn);

% 1. set number of folds
%[c b] = hist( labels,1:1:max(labels));
[c b] = hist( labels,unique(labels));
idx = find(c==0);
c(idx) = [];
b(idx) = [];
mi = min(c);

if ~exist( 'N','var')
    N = 10;
end
if mi<N
    N=mi;
end

% 2. split dataset into N folds
splits = [];
splits{N} = [];
for i=1:length(b)
    idx = find( labels==b(i));
    rp = randperm(length(idx));  % randomize partitioning
    val = floor(length(idx)/N);
    array = val+zeros(N,1);
    re = length(idx)-val*N;
    array(1:re) = array(1:re)+1;
    st = 1;
    for j=1:N
        en = st+array(j);
        splits{j} = [splits{j}; idx(rp(st:en-1))];
        st = en;
    end
end
