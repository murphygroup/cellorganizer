function yAll = compare_shape_model_spharm_obj(models, param, fileID)

idx = [0];
xAll = [];
for i=1:length(models)
    xAll = [xAll;models{i}.X];
    idx = [idx idx(end)+size(models{i}.X, 1)];
end
[scales,coeff,score,latent,tsquared,explained,mu,train_score_consolidate,train_explained,yAll] = CalculatePCA(xAll, 3);
yAll = [];
for i=1:length(models)
    yAll{i} = train_score_consolidate(idx(i)+1:idx(i+1), :);
end
    
for j=1:length(models)-1
    pool = vertcat(yAll{j}, yAll{end});
    d1 = idx(j+1) - idx(j);
    d2 = idx(end) - idx(end-1);
    prob2 = d1/(d1+d2);
    temp = zeros(1, d1+d2);
    for i=1:(d1+d2)
        Idx = knnsearch(pool,pool(i, :), 'K', 9);
        prob1 = sum(Idx(2:end) < d1)/8;
        temp(i)=prob1.*log((prob1 + 1e-7)/prob2) + (1-prob1).*log((1-prob1+ 1e-7)/(1-prob2));
    end


    text2html(fileID, sprintf('Shape space comparison score between model1 and model2: %.4f;\n', mean(temp)));
    p_val = permutation_test(mean(temp), yAll{end});
    text2html(fileID, sprintf('P value of shape space comparison score: %.4f;\n', p_val));
    % KLD = clique_percolation_spatial(y1, y2);
    % text2html(fileID, sprintf('Clique percolation divergence between model1 and model2: %.4f;\n', KLD));
  

end
end



