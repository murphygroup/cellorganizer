function [d_source, d_target, d_source_distance, d_target_distance, r_velocity_offset] = ...
  Greedy3D_single_step_window(windowed_source, windowed_target, velocity_offset, options)
  % Internal function for Greedy3D_lambda_pre_compressed.
  %
  % Computes the derivatives of the deformation fields and distances for particular windows of the (currently deformed) source and target images.
  % 
  % 2013-02-17 tebuck: Changing the default options to those that are currently most useful or are being used in CellOrganizer.


  % Default option values and options processing (these are now set in Greedy3D_compressed_default_options.m):

  % default_options.alpha = 1;          % coefficient of Lapacian operator
  % default_options.gamma = 0.04;       % Coefficient of the Identity
  % default_options.single_sided = false;
  % default_options.use_fft = true;
  % default_options.use_gaussian_kernel = false;
  % default_options.drop_kernel = true;
  % default_options.kernel_z_radius = 32;
  % default_options.gradient_type = 'scharr5';
  % default_options.scale_velocities_like_kernel = true;
  % default_options.degenerate_kernel = false;
  % default_options.gradient_mitering = false;
  % default_options.convolution_method = 'convnfft_fast';
  % default_options.filter_radius = 32; 
  default_options = Greedy3D_compressed_default_options();
  if ~exist('options', 'var')
    options = default_options;
  else
    options = process_options_structure(default_options, options);
  end
  
  have_velocity_offset = exist('velocity_offset', 'var');
  if (have_velocity_offset)
    have_velocity_offset = ~isempty(velocity_offset); 
  end
  
  

  S = windowed_source;
  src2 = S;
  T = windowed_target;
  trg2 = T;
  source_size = size(S);
  source_size(numel(source_size) + 1:3) = 1; 
  %filter_size = size(S) .* [0, 0, 1] + [2, 2, 0] * options.filter_radius;
  filter_size = options.kernel_z_radius .* [0, 0, 2] + [2, 2, 0] * options.filter_radius;
  M = source_size(1); 
  N = source_size(2); 
  P = source_size(3); 
  M2 = filter_size(1); 
  N2 = filter_size(2); 
  P2 = filter_size(3); 

  valid_inds_x = (options.filter_radius + 2):(N - options.filter_radius - 1);
  valid_inds_y = (options.filter_radius + 2):(M - options.filter_radius - 1);
  valid_inds_z = 1:P;
  
  is2d = P == 1; 
  
  src2(isnan(src2))=0; 
  if (~options.single_sided)
    trg2(isnan(trg2))=0;  
  end
  
  kernel = [];
  if ~options.degenerate_kernel
    if ~options.use_gaussian_kernel
      kernel = Greedy3D_kernel_resampled([M2, N2, P2], options.alpha, options.gamma, true);
      kernel = real(kernel);
      if options.drop_kernel
        kernel_sum = abs(sum(kernel(:)));
        kernel = kernel - max(kernel(:));
        kernel = kernel .* (kernel_sum ./ abs(sum(kernel(:))));
      end
    else
      kernel = -Greedy3D_kernel_gaussian([M2, N2, P2], options.alpha, options.gamma);
    end
  end
  
  %compute force F = -(T - src).*grad(src)
  [Fx1,Fy1,Fz1] = smooth_gradient(src2, options.gradient_type, options.gradient_mitering);
  [Fx2,Fy2,Fz2] = smooth_gradient(trg2, options.gradient_type, options.gradient_mitering);
  im_diff = trg2 - src2;% Computing the difference in images
  clear src2 trg2
  Fx1 = -im_diff.*Fx1;
  Fy1 = -im_diff.*Fy1;		
  Fz1 = -im_diff.*Fz1;		
  Fx2 = im_diff.*Fx2;
  Fy2 = im_diff.*Fy2;
  Fz2 = im_diff.*Fz2;
	
  % Contribution of the source image
  % Computing the velocity by solving (LtL)v = F for v
  if ~options.degenerate_kernel
    if ~options.use_fft
      Vx1 = imfilter(Fx1, kernel, 'replicate', 'same', 'conv');
      Vy1 = imfilter(Fy1, kernel, 'replicate', 'same', 'conv');
      Vz1 = imfilter(Fz1, kernel, 'replicate', 'same', 'conv');
    else
      switch lower(options.convolution_method)
        case 'convnfft'
%           options.Power2Flag = true;
          options.GPU = false;
          Vx1 = convnfft(Fx1, kernel, 'same', [], options);
          Vy1 = convnfft(Fy1, kernel, 'same', [], options);
          Vz1 = convnfft(Fz1, kernel, 'same', [], options);
        case 'convn'
          Vx1 = convn(Fx1, kernel, 'same');
          Vy1 = convn(Fy1, kernel, 'same');
          Vz1 = convn(Fz1, kernel, 'same');
        case 'convolution3d_fftdomain'
          Vx1 = convolution3D_FFTdomain(Fx1, kernel);
          Vy1 = convolution3D_FFTdomain(Fy1, kernel);
          Vz1 = convolution3D_FFTdomain(Fz1, kernel);
        otherwise
          Vx1 = convnfft_fast(Fx1, kernel);
          Vy1 = convnfft_fast(Fy1, kernel);
          Vz1 = convnfft_fast(Fz1, kernel);
      end
    end
  else
    Vx1 = -Fx1;
    Vy1 = -Fy1;
    Vz1 = -Fz1;
  end
  
  if (~have_velocity_offset)
    velocity_offset = [...
      Vx1(valid_inds_y(1), valid_inds_x(1), valid_inds_z(1))...
      Vy1(valid_inds_y(1), valid_inds_x(1), valid_inds_z(1))...
      Vz1(valid_inds_y(1), valid_inds_x(1), valid_inds_z(1))...
                      ];
  end
  r_velocity_offset = velocity_offset;
  Vx1 = Vx1 - velocity_offset(1);  % Specify arbitrary constant by forcing corner u = 0.
  Vy1 = Vy1 - velocity_offset(2); 
  Vz1 = Vz1 - velocity_offset(3); 

  Vx1(:,1,:) = 0;Vx1(:,N,:) = 0;Vx1(1,:,:) = 0;Vx1(M,:,:) = 0;
  Vy1(:,1,:) = 0;Vy1(:,N,:) = 0;Vy1(1,:,:) = 0;Vy1(M,:,:) = 0;
  Vz1(:,1,:) = 0;Vz1(:,N,:) = 0;Vz1(1,:,:) = 0;Vz1(M,:,:) = 0;
  if ~is2d
      Vx1(:,:,1)=0;Vx1(:,:,P)=0;
      Vy1(:,:,1)=0;Vy1(:,:,P)=0;
      Vz1(:,:,1)=0;Vz1(:,:,P)=0;
  else
    Vz1 = zeros(size(Vz1));
  end
  
  Vx1(isnan(Vx1))=0; 
  Vy1(isnan(Vy1))=0; 
  Vz1(isnan(Vz1))=0; 

  clear Fx1 Fy1 Fz1 fhat_x fhat_y fhat_z
  
  if options.scale_velocities_like_kernel
    relative_scales = filter_size ./ 64.; 
    Vx1 = Vx1 * relative_scales(2);
    Vy1 = Vy1 * relative_scales(1);
    Vz1 = Vz1 * relative_scales(3);
  end
  
  % Computing the distance by (1/N*M)*sqrt(dx2+ dy2);
  dx2 = (-del2(Vx1)*6*options.alpha + options.gamma*Vx1).^2;
  dy2 = (-del2(Vy1)*6*options.alpha + options.gamma*Vy1).^2;
  dz2 = (-del2(Vz1)*6*options.alpha + options.gamma*Vz1).^2;
  
  dx2 = mean2(dx2(valid_inds_y, valid_inds_x, valid_inds_z));
  dy2 = mean2(dy2(valid_inds_y, valid_inds_x, valid_inds_z));
  dz2 = mean2(dz2(valid_inds_y, valid_inds_x, valid_inds_z));
  ddist1 = sqrt(dx2 + dy2 + dz2);
  Vx1 = Vx1(valid_inds_y, valid_inds_x, valid_inds_z);
  Vy1 = Vy1(valid_inds_y, valid_inds_x, valid_inds_z);
  Vz1 = Vz1(valid_inds_y, valid_inds_x, valid_inds_z);

  
  d_source = {Vx1, Vy1, Vz1};
  d_source_distance = ddist1;
  if (options.single_sided)
    % We're done:
    d_target = [];
    d_target_distance = [];
    return
  end
  
  
  % contribution of target image
  
  if ~options.degenerate_kernel
    if ~options.use_fft
      Vx2 = real(imfilter(Fx2, kernel, 'replicate', 'same', 'conv'));
      Vy2 = real(imfilter(Fy2, kernel, 'replicate', 'same', 'conv'));
      Vz2 = real(imfilter(Fz2, kernel, 'replicate', 'same', 'conv'));
    else
      switch lower(options.convolution_method)
        case 'convnfft'
          options.Power2Flag = false;
          Vx2 = convnfft(Fx2, kernel, 'same', [], options);
          Vy2 = convnfft(Fy2, kernel, 'same', [], options);
          Vz2 = convnfft(Fz2, kernel, 'same', [], options);
        case 'convn'
          Vx2 = convn(Fx2, kernel, 'same');
          Vy2 = convn(Fy2, kernel, 'same');
          Vz2 = convn(Fz2, kernel, 'same');
        case 'convolution3d_fftdomain'
          Vx2 = convolution3D_FFTdomain(Fx2, kernel);
          Vy2 = convolution3D_FFTdomain(Fy2, kernel);
          Vz2 = convolution3D_FFTdomain(Fz2, kernel);    
              
        otherwise
          Vx2 = convnfft_fast(Fx2, kernel);
          Vy2 = convnfft_fast(Fy2, kernel);
          Vz2 = convnfft_fast(Fz2, kernel);
      end
    end
  else
    Vx2 = -Fx2;
    Vy2 = -Fy2;
    Vz2 = -Fz2;
  end
  
  if (~have_velocity_offset)
    velocity_offset = [...
      velocity_offset ...
      Vx2(valid_inds_y(1), valid_inds_x(1), valid_inds_z(1))...
      Vy2(valid_inds_y(1), valid_inds_x(1), valid_inds_z(1))...
      Vz2(valid_inds_y(1), valid_inds_x(1), valid_inds_z(1))...
                      ];
  end
  r_velocity_offset = velocity_offset;
  Vx2 = Vx2 - velocity_offset(4);  % Specify arbitrary constant by forcing corner u = 0.
  Vy2 = Vy2 - velocity_offset(5); 
  Vz2 = Vz2 - velocity_offset(6); 

  Vx2(:,1,:) = 0;Vx2(:,N,:) = 0;Vx2(1,:,:) = 0;Vx2(M,:,:) = 0;
  Vy2(:,1,:) = 0;Vy2(:,N,:) = 0;Vy2(1,:,:) = 0;Vy2(M,:,:) = 0;
  Vz2(:,1,:) = 0;Vz2(:,N,:) = 0;Vz2(1,:,:) = 0;Vz2(M,:,:) = 0;
  if ~is2d
    Vx2(:,:,1)=0;Vx2(:,:,P)=0;
    Vy2(:,:,1)=0;Vy2(:,:,P)=0;
    Vz2(:,:,1)=0;Vz2(:,:,P)=0;
  else
    Vz2 = zeros(size(Vz2));
  end
  
  Vx2(isnan(Vx2))=0; 
  Vy2(isnan(Vy2))=0; 
  Vz2(isnan(Vz2))=0; 

  clear Fx2 Fy2 Fz2 fhat_x fhat_y fhat_z
  
  if options.scale_velocities_like_kernel
    relative_scales = filter_size ./ 64.; 
    Vx2 = Vx2 * relative_scales(2);
    Vy2 = Vy2 * relative_scales(1);
    Vz2 = Vz2 * relative_scales(3);
  end

  dx22 = (-del2(Vx2)*6*options.alpha + options.gamma*Vx2).^2;
  dy22 = (-del2(Vy2)*6*options.alpha + options.gamma*Vy2).^2;
  dz22 = (-del2(Vz2)*6*options.alpha + options.gamma*Vz2).^2;

  dx22 = mean2(dx22(valid_inds_y, valid_inds_x, valid_inds_z));
  dy22 = mean2(dy22(valid_inds_y, valid_inds_x, valid_inds_z));
  dz22 = mean2(dz22(valid_inds_y, valid_inds_x, valid_inds_z));
  ddist2 = sqrt(dx22 + dy22 + dz22);
  Vx2 = Vx2(valid_inds_y, valid_inds_x, valid_inds_z);
  Vy2 = Vy2(valid_inds_y, valid_inds_x, valid_inds_z);
  Vz2 = Vz2(valid_inds_y, valid_inds_x, valid_inds_z);
  
  clear dx2 dy2 dz2 dx22 dy22 dz22
  
  d_target = {Vx2, Vy2, Vz2};
  d_target_distance = ddist2;
  
  clear Vx1 Vy1 Vz1 Vx2 Vy2 Vz2
  

