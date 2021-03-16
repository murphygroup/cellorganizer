function all_script(cellnum,batchno)

rand_start = 1;
%subfolder = 1;  % Aabid's Gaussian model
%subfolder = 2;  %%Jieyue, Erlang, one-step rebound control
subfolder = 3;  %%Jieyue, Erlang, multi-steps rebound control

%n_inputArray = [5,50:50:350];
%mulen_inputArray  = [5 10 15 20 25 30 35];
n_inputArray = [5, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 350];
mulen_inputArray  = [5 10 15 20 25 30 35 40];
colli_inputArray = [0.9 0.95 0.98];

if subfolder==1
coeff_var_Array = [0 0.1 0.2 0.3];
else
coeff_var_Array = [0];
end

extract_features_hpa_for_part2_red(subfolder,cellnum,batchno,n_inputArray,mulen_inputArray,colli_inputArray,coeff_var_Array);
