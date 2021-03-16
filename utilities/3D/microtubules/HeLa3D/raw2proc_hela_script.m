clear all
close all

% some cells in here were left out because they have imaging errors

for I = [1:4 6:42 44:52]
	raw2proc_hela(I);
end
