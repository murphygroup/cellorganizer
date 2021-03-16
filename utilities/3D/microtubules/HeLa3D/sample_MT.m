function [n,mu_len,collin] = sample_MT(resultpath)
%created by Devin Sullivan 2/18/12
%
%
%This code samples from the values of the HeLa dataset with replacement for
%the three parameters (number, mu_len, collin)
%
%Inputs:
%resultpath = path to mat file containing results for the dataset
load('HeLaresult.mat')