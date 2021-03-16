function [wopt,wtot] = classweights( classlabels)

wopt = ' ';
[U I J] = unique(classlabels);
[stot b] = hist(J,1:1:max(J));
wtot = min(stot)./stot;
for j=1:length(wtot)
    wopt = [wopt '-w' num2str(U(j)) ' ' num2str(wtot(j)) ' '];
end

