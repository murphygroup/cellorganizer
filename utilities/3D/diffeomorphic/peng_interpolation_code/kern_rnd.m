function xx = kern_rnd(x,h)
%KERN_RND draw random number from a kernel distribution with Gaussian
%kernel of bandwidth H

[n,d] = size(x);
Sigma = (h^2)*eye(d);
randseq = randperm(n);
i = randseq(1);
mu = x(i,:);
xx = mvnrnd(mu,Sigma);