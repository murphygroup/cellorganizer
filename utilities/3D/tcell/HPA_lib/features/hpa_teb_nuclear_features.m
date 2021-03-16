function [features, names, slf_names] = ...
    hpa_teb_nuclear_features(dna_image, nuclear_region)
% hpa_teb_nuclear_features(dna_image, nuclear_region)
%%    hpa_teb_nuclear_features(dna_image, nuclear_region, cellular_region)
% Created: 2009-03-26 20:00 tebuck 
% Assume background-subtracted! No more bg sub prior to feat calc.

  % Estimated by Justin Newberg, microns per pixel side:
  latres = .05;
  n = {}; 
  sn = {}; 
  number_haralick_resolutions = 5; 
  %number_haralick_resolutions = 4; 
  f = ones(0, 0);

    % Haralick features on the nucleus
    [n2, f2, sn2] = ml_featset( ...
      dna_image, ...
      nuclear_region, [], 'har', ...
      latres, [], 'nobgsub', 'none' , latres, 32 ...
      ); 
    %  latres, [], 'lowcommon', 'nih' , latres, 32 ...
    n2 = strcat(n2, [' (nuclear)']); 
    sn2 = strcat(sn2, [' (nuclear)']); 
    n = [n, n2]; 
    sn = [sn, sn2]; 
    f = [f, f2]; 
    
    for j = 1: (number_haralick_resolutions - 1)
      % Haralick features on the nucleus
      [n2, f2, sn2] = ml_featset( ...
        dna_image, nuclear_region, [], 'har', ...
        latres, [], 'nobgsub', 'none' , latres*(2^j), 32 ...
        ); 
      n2 = strcat(n2, [' (nuclear, res' num2str(1 / (2^j)) ')']); 
      sn2 = strcat(sn2, [' (nuclear, res' num2str(1 / (2^j)) ')']); 
      n = [n, n2]; 
      sn = [sn, sn2]; 
      f = [f, f2]; 
    end
    
    [f2, n2, sn2] = hpa_psdnafeats( ...
      dna_image, ...
      nuclear_region, latres); 
    %  ml_imgbgsub(dna_image, 'lowcommon'), ...
    %  nuclear_region, cellular_region, latres); 
    n = [n, n2]; 
    sn = [sn, sn2]; 
    f = [f, f2]; 
    
    features = f; 
    names = n; 
    slf_names = sn; 
    

    
function [ r_features, r_feature_names, r_feature_slfnames ] = ...
    hpa_psdnafeats ( p_nuclear_image, p_nuclear_crop, ...
                                     resolution)

% Created: 2009-03-27 14:52 TEB, copied from get_cell_psdnafeats.m
% Operates on single cell with crop image
% Assume background-subtracted! No more bg sub prior to feat calc.

r_features = [];
r_feature_names = cell(1, 0);
r_feature_slfnames = cell(1, 0);

    nucl_img = uint16(double(p_nuclear_image) .* double(p_nuclear_crop));

        nucl_thresh = p_nuclear_crop;
        nucl_proc = nucl_img;
        nucl_area = sum(sum(nucl_thresh)) * (resolution * resolution);

    temp_features = zeros(1, 5);
            r_feature_names{1} = 'total_nuclear_intensity';
            r_feature_names{2} = 'total_nuclear_area (square micrometers)';
            r_feature_names{3} = 'nuclear_intensity_times_area (times square micrometers)';
            r_feature_names{4} = 'mean_nuclear_intensity (per square micrometer)';
            r_feature_names{5} = 'nuclear_compactness';
      r_feature_slfnames = repmat({'-'}, 1, 5);

            % Calculate total DNA, total area above zero,
            % and DNA per pixel:
            % The prot total will be over the region and so

            temp_features(1) = sum(sum(nucl_proc));
            temp_features(2) = nucl_area;
            temp_features(3) = temp_features(1) * temp_features(2);
            temp_features(4) = 0;
            if ( temp_features(2) > 0 )
                temp_features(4) = temp_features(1) / temp_features(2);
            end;

            nucl_mask_er = logical(nucl_thresh);
            nucl_comp = obj_compactness(nucl_mask_er);
            temp_features(5) = nucl_comp;

        r_features = temp_features;

    

        
function [ r_C ] = obj_compactness( p_img )
P = obj_perimeter4(p_img);
A = obj_area(p_img);
if ( P * A == 0 )
	r_C = 0;
	return;
end;
r_C = (4 * pi * A) / P^2;

function [ r_A ] = obj_area( p_img )
uvc = length(unique(p_img(:)));
if ( uvc <= 1 ); r_A = 0; return; end;
if ( uvc == 2 )
	p_img = double(p_img);
end;
[A] = regionprops(p_img, 'Area');
if ( isempty(A) )
	r_A = 0;
	return;
end;
r_A = zeros(1, length(A));
for i=1:length(A)
	r_A(i) = A(i).Area;
end;

function [ r_P ] = obj_perimeter4( p_img )
  
  width = size(p_img, 2); 
  height = size(p_img, 1); 
  p_img = [zeros(1, width + 2); ...
           zeros(height, 1), p_img, zeros(height, 1); ...
           zeros(1, width + 2)]; 
  r_P = sum(sum(bwperim(p_img, 4)));
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
