function tp_epcdf(samples)
% Plot the experimental CDF of a group of data

x = sort(unique(samples));
H = hist(samples,x);
plot(x,cumsum(H)/length(samples),'k.','MarkerSize',2.5);