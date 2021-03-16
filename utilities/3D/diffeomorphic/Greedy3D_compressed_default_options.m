function [default_options] = Greedy3D_compressed_default_options()
  % Internal function for Greedy3D_lambda_pre_compressed.m.
  % 
  % Provides a set of options that are the defaults across Greedy3D_lambda_pre_compressed.m and its internal functions.
  %
  % Currently most useful options:
  % known_distance, stop_source_distance (useful if interpolating and the (approximate) true distance between the input shapes and between the source and the desired shape are known)
  % single_sided (useful if registering one image to another)
  % just_compute_distance (useful if not interested in the interpolated image)
  % maximum_deformation_per_step (useful for images of thin shapes)
  % window_radius (always set this to the width of your images unless you want windowing)
  % filter_radius, kernel_z_radius,scale_velocities_like_kernel (change these to match the aspect ratio of and size of interesting variations in your images)
  % 
  default_options = struct();
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Used in Greedy3D_lambda_pre_compressed:

  % Default option values and options processing:
  % known_distance says how far the shapes are known/assumed to be (negative values disable this function):
  default_options.known_distance = -1;
  % stop_source_distance is how far the source can move before stopping, and (known_distance - stop_source_distance) is how far the target can move before stopping (negative values disable this function):
  default_options.stop_source_distance = -1;
  % If true, registering one image to another rather than symmetrically interpolating them:
  default_options.single_sided = false; 
  % Size of the regularization filter/kernel used each iteration of LDDMM:
  default_options.filter_radius = 32;
  % Size of windows (actually the full diameter):
  default_options.window_radius = 64;
  % If true, windows will be compressed in memory using GZIP:
  default_options.use_compression = false; 
  % If true, no interpolation result is returned because just the distance between images is needed as output:
  default_options.just_compute_distance = false; 
  % If true, save a .mat file for each iteration's result:
  default_options.save_intermediates = false; 
  % Specifies that Greedy3D_lambda_pre_compressed.m's lambda_list contains distances rather than distance ratios (see Greedy3D_lambda_pre_compressed.m's description of lambda_list):
  default_options.lambdas_are_distances_from_source = false; 

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Used in Greedy3D_integration_compressed:
  
  % Runge-Kutta integrator parameters (Bogacki-Shampine method):
  % Step size fractions:
  default_options.integrator_c = [0, 1/2, 3/4, 1];
  % Stage coefficients for the two final estimates (b is the integration formula, b2 the error/comparison formula):
  default_options.integrator_b = [2/9, 1/3, 4/9, 0];
  default_options.integrator_b2 = [7/24, 1/4, 1/3, 1/8];
  % Per-stage coefficients for the linear combination of prior stages:
  default_options.integrator_a = zeros(4, 3);
  default_options.integrator_a(2, 1) = 1/2;
  default_options.integrator_a(3, 1:2) = [0, 3/4];
  default_options.integrator_a(4, 1:3) = [2/9, 1/3, 4/9];
  % Tells us whether the last stage of one iteration is the same as
  % the first stage of the next iteration:
  default_options.integrator_first_same_as_last = true;
  % Order of error in integration results (is this right?):
  default_options.integrator_error_order = 3;
  
  % Initial integration step size:
  default_options.step_size = 0.01;
  
  % Convergence criteria:
  % Threshold of distance moved per unit time below which to declare convergence:
  default_options.convergence_tolerance = 0;
  % Threshold of sum of squared errors between the deformed source and target pixels:
  default_options.convergence_registration_error = 0;
  % Threshold of sum of absolute errors between the deformed source and target pixels:
  default_options.convergence_absolute_error = 0;
  % Threshold of change in sum of absolute errors per iteration below which to declare convergence:
  default_options.convergence_absolute_error_difference = 0;
  % Number of iterations over which to measure absolute errors for iteration_convergence_absolute_error:
  default_options.convergence_absolute_error_difference_window_size = 5;
  
  % % known_distance says how far the shapes are known/assumed to be (negative values disable this function):
  % default_options.known_distance = -1;
  % % stop_source_distance is how far the source can move before stopping, and (known_distance - stop_source_distance) is how far the target can move before stopping (negative values disable this function):
  % default_options.stop_source_distance = -1;
  
  % Step size control criteria:
  % Maximum LDDMM distance error (between the two integrations, if applicable) change allowed per iteration:
  default_options.absolute_distance_tolerance = 0; 
  % Maximum deformation error in any coordinate for any pixel allowed per iteration:
  default_options.absolute_deformation_tolerance = 0; 

  % Switches for saving intermediate results:
  % If true, save a .mat file for each iteration's result:
  % default_options.save_intermediates = false;
  % If true, save an image showing the current source and target images for each iteration's result:
  default_options.save_intermediate_images = false;
  % If true, save the source, target, and deformed source and target when save_intermediates is true:
  default_options.always_save_full_intermediates = false; 
  % If true, save each RK stage's result to disk, otherwise store the last two iterations in memory (for interpolation, requires a second integration unless known_distance and stop_source_distance are provided):
  default_options.save_stages = false;

  % % No interpolation need be returned if just computing distance between images:
  % default_options.just_compute_distance = false; 
  % % True if registration, false if interpolation:
  % default_options.single_sided = false; 
  % Convergence should not be dependent on number of iterations, so set this to a large number:
  default_options.maximum_iterations = 10000; 
  % If the image space is not toroidal/periodic, don't use periodic convolution so, e.g., the left and right sides of images do not interact (if this is false):
  default_options.periodic_space = false; 
  % Old option that quits if the registration error goes up too many times (disabled if 0):
  default_options.maximum_registration_error_failures = 0;
  % Old option that multiplies the step size by this amount, if greater than 0, upon registration error failure:
  default_options.registration_error_failure_step_scale = 0;
  % The amount of spatial deformation along each coordinate permitted per iteration. Useful if cells are very flat, in which case large movements in the normal direction during integration might fold the cell over itself (i.e., not diffeomorphic):
  default_options.maximum_deformation_per_step = ones(1, 3) * .5; 
  % If true, print statistics during integration:
  default_options.verbose = 1; 
  % 1x2 double array, the time over which to integrate the differential equations that deform the input images if convergence conditions are not met (set to [0, inf] by default because now there are convergence criteria; by the way, max_time has been ignored for a long time):
  default_options.time_span = [0, inf]; 

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Used in Greedy3D_single_step_compressed:
  
  % default_options.single_sided = false; 
  % default_options.periodic_space = false; 
  % default_options.filter_radius = 32; 
  % default_options.verbose = 1; 
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Used in Greedy3D_single_step_window:
  
  default_options.alpha = 1;          % coefficient of Lapacian operator
  default_options.gamma = 0.04;       % Coefficient of the Identity
  % default_options.single_sided = false;
  % If true, do convolutions in the frequency domain for efficiency:
  default_options.use_fft = true;
  % If true, use a Gaussian as the low pass filter instead of the squared inverse of the differential operator paper from Rohde et al. 2008 (i.e., (alpha * del^2 + gamma)^-2):
  default_options.use_gaussian_kernel = false;
  % If true, approximate the apparent limit of the Laplacian-based low pass filter as image size goes to infinity by setting the minimum value to zero:
  default_options.drop_kernel = true;
  % Like filter_radius, but in the third dimension:
  default_options.kernel_z_radius = 32;
  % See smooth_gradient for how to choose gradient_type:
  default_options.gradient_type = 'scharr5';
  % If true, velocities (i.e., deformations at each iteration) will be larger in the same directions for which the kernel is larger:
  default_options.scale_velocities_like_kernel = true;
  % If true, do not use a low pass filter/regularize deformations:
  default_options.degenerate_kernel = false;
  % If true, attempt to correct the cancellation of gradients in sharp corners (this did not seem to prevent corners from lagging behind smoother edges):
  default_options.gradient_mitering = false;
  % If convnfft_fast, use our faster version of convnfft that excludes many of the options from the File Exchange version of that function:
  default_options.convolution_method = 'convnfft_fast';
  % default_options.filter_radius = 32; 
  default_options.tempparent = 'temp';
  
  
  
  
  
  
