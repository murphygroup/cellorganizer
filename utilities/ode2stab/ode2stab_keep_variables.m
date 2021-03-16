function [some_t, some_x, all_t, all_x] = ...
    ode2stab_keep_variables(...
      odesolver,odefunc,tspan,Y0,tchunk,maxcv...
      , keep_proportion, keep_variables, options...
      )
% ODE2STAB_KEEP_INDICES Integrates systems of differential equations until stability is
% reached (or time is over), including arguments to discard some
% integration variables for memory conservation
% Original help:
%	[TOUT,YOUT] = ODE2STAB(ODESOLVER,ODEFUNC,TSPAN,Y0,TCHUNK,MAXCV) with 
%	TSPAN = [T0 TFINAL], integrates the system of differential equations 
%	stored in file ODEFUNC until stability is reached or until the time 
%	specified in TFINAL, whatever happens first. Stability is specified 
%	as a maximum coefficient of variation MAXCV for all components of the 
%	solution. The first argument (odesolver) must be a function handle for one of the
%	MATLAB ODE solvers, i.e. @ode45. Y0 is a vector with the initial
%	conditions.
%	The function integrates the system in 'chunks' of time specified in 
%	TCHUNK, which should be sufficiently smaller than TFINAL-T0. After each 
%	integration run, the two conditions are chequed. If the coefficient of
%	variation of all the components of the solution is smaller than MAXCV
%	OR we have reached TFINAL, the integration stops.
%	Example:
%
% 	First create a file called odelogist.m as:
%	----------------------------------
% 	function dn= odelogist(t,x)
% 	dn= zeros(1,1);
% 	dn(1)= (0.01*(1-(x(1)/10)))*x(1);
%	----------------------------------
%
% 	odefunc= 'odelogist';
% 	tspan= 1:10000;
% 	Y0= 0.01;
% 	tchunk= 100;
% 	maxcv= 10^-3;
% 	[t,x]= ode2stab(@ode45,odefunc,tspan,Y0,tchunk,maxcv);
% 	plot(t,x)
% 
%	Integration stops at time 1101, instead of 10000.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Original code under BSD license (see ode2stab's license as of 2011-05-19).
%
% Modified 2011-05-12 15:31 by tebuck to optionally discard intermediate results.
% Modified 2011-05-13 12:46 by tebuck to optionally keep the values
% of some variables for all time points (as well as the times).
% Modified 2011-05-19 by tebuck to optionally prevent the solver
% from keeping intermediate values by using the light_integration
% option in the options structure. It would be nice if the
% interface were more consistent, i.e., this option were given as
% an argument, but it's this way for now.
% 2011-05-23 16:04 tebuck: This version should return all the time
% points in tspan in addition to its previous functionality.

%tspan
tspan = unique(tspan); 

some_x= [];
some_t= [];
all_x= [];
all_t= [];
cv= 2*maxcv;
totaltime= 0;
%time= (tspan(1):tspan(1)+tchunk)';
%time= [tspan(1), tspan(1)+tchunk];
%time= [tspan(1), min(tspan(1) + tchunk, tspan(end))];
next_time_index = 2; 
next_time = tspan(next_time_index); 
time= [tspan(1), min(tspan(1) + tchunk, next_time)];
x= Y0;
end_of_time = false; 

%-- Check that we have a function handle for the solver
if ~isa(odesolver,'function_handle') 
	error('First argument must be a function handle of a ODE solver, i.e. @ode45');
	return
end


%-- Run the thing until maxtime or until CV is less than requested
stop_condition = false; 
%while cv > maxcv & totaltime < tspan(length(tspan))
start_wall_time = tic;
while ~stop_condition
	ninit= x(size(x,1),:);
	%[t,x]= odesolver(odefunc,time,ninit);
  %if (options.light_integration)
  if (exist('options.light_integration', 'var') && options.light_integration)
    time = [time(1), mean(time), time(2)]; 
  end
  time
	[t,x]= odesolver(odefunc,time,ninit, options);
	totaltime= totaltime + time(end) - time(1);
	%time= (time(length(time))+1:time(length(time))+tchunk)';
	%time= [time(2), time(2)+tchunk];
	%time= [time(2), min(time(2) + tchunk, tspan(end))];
  %if (time(end) == next_time && next_time_index < length(tspan))
  if time(end) >= next_time
    if next_time_index < length(tspan)
      next_time_index = next_time_index + 1; 
      next_time = tspan(next_time_index); 
    else
      end_of_time = true; 
    end
  end
  
	time= [time(end), min(time(end) + tchunk, next_time)];
	%totaltime= totaltime+ tchunk;
	cv= max(var(x)./mean(x));
  if (length(keep_variables) > 0)
    %all_x = [all_x; x(:, keep_variables)];
    %all_t = [all_t; t];
    all_x = [all_x; x(1:(end - 1), keep_variables)];
    all_t = [all_t; t(1:(end - 1), :)];
  end
  keep_indices = rand(size(x, 1), 1) <= keep_proportion;
  if (keep_proportion == 0)
    keep_indices = false(size(keep_indices)); 
  end
  %arrayfun(@(y)(max(y == tspan)), t)
  %keep_indices
  keep_indices = keep_indices | arrayfun(@(y)(max(y == tspan)), t);
  keep_indices(end) = false;
  
  if (sum(keep_indices) > 0)
    some_x= [some_x; x(keep_indices, :)];
    some_t= [some_t; t(keep_indices, :)];
  end
  %stop_condition = ~(cv > maxcv & totaltime < tspan(length(tspan)));
  stop_condition = ~(cv > maxcv && totaltime < tspan(length(tspan)) ...
                     && ~end_of_time);
  if stop_condition
% $$$     if keep_proportion == 0
% $$$       some_x = x(end, :); 
% $$$       some_t = t(end, :); 
% $$$     else
% $$$       some_x = [some_x; x(end, :)]; 
% $$$       some_t = [some_t; t(end, :)]; 
% $$$     end
    some_x = [some_x; x(end, :)]; 
    some_t = [some_t; t(end, :)]; 
    all_x = [all_x; x(end, keep_variables)];
    all_t = [all_t; t(end, :)];
  end
	clear t
  %toc
  fprintf('%f time elapsed in %s\n', toc(start_wall_time), mfilename);
end
%toc