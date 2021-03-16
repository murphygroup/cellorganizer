function [psf] = which_psf(microscope) 
%Devin Sullivan 2/18/12
%
%
%This function takes in a string name of a microscope for which to return a
%psf for, if 'none' is specified this will return an empty matrix 
%
%Input:
%microscope = string name for the microscope 
%
%
%Output:
%psf = name of .jpg file for the psf for the specified microscope contained
%in the HeLa3D/psf folder

switch microscope
    case 'SVI'
        psf = 'HPA_285.mat';
%     case 'Andor'
%         psf = 'example.jpg';
    case 'none'
        psf = [];
    otherwise
        disp('Error, no microscope specified, applying no psf!');
        psf = [];
end