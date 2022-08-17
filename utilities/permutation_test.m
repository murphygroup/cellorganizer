function p_val = permutation_test(obs_score, X1)
[d,~] = size(X1);

d1 = d/2;
d2 = d - d1;


score_list = [];
for k=1:50
    idx1 = randperm(d);
    X = X1(idx1, :);

    prob2 = d1/(d1+d2);
    temp = zeros(1, d1+d2);
    for i=1:(d1+d2)
        Idx = knnsearch(X, X(i, :), 'K', 9);
        prob1 = sum(Idx(2:end) < d1)/8;
        temp(i)=prob1.*log((prob1 + 1e-7)/prob2) + (1-prob1).*log((1-prob1+ 1e-7)/(1-prob2));
    end
    score_list(end+1) = mean(temp);
end
% load golgi_list.mat score_list
p_val = sum(obs_score<score_list)/50;
end