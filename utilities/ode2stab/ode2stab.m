function [allt,allx]= ode2stab(odesolver,odefile,tspan,Y0,tchunk,maxcv...
                               , keep_proportion, options...
                               )
% ODE2STAB Integrates systems of differential equations until stability is
% reached (or time is over)
%	[TOUT,YOUT] = ODE2STAB(ODESOLVER,ODEFILE,TSPAN,Y0,TCHUNK,MAXCV) with 
%	TSPAN = [T0 TFINAL], integrates the system of differential equations 
%	stored in file ODEFILE until stability is reached or until the time 
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
% 	odefile= 'odelogist';
% 	tspan= 1:10000;
% 	Y0= 0.01;
% 	tchunk= 100;
% 	maxcv= 10^-3;
% 	[t,x]= ode2stab(@ode45,odefile,tspan,Y0,tchunk,maxcv);
% 	plot(t,x)
% 
%	Integration stops at time 1101, instead of 10000.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Original code under BSD license.
%
% Modified 2011-05-12 15:31 by tebuck to optionally discard intermediate results.


allx= [];
allt= [];
cv= 2*maxcv;
totaltime= 0;
%time= (tspan(1):tspan(1)+tchunk)';
time= [tspan(1), tspan(1)+tchunk];
x= Y0;

%-- Check that we have a function handle for the solver
if ~isa(odesolver,'function_handle') 
	error('First argument must be a function handle of a ODE solver, i.e. @ode45');
	return
end


%-- Run the thing until maxtime or until CV is less than requested
stop_condition = false; 
%while cv > maxcv & totaltime < tspan(length(tspan))
tic
while ~stop_condition
	ninit= x(size(x,1),:);
	%[t,x]= odesolver(odefile,time,ninit);
	[t,x]= odesolver(odefile,time,ninit, options);
	%time= (time(length(time))+1:time(length(time))+tchunk)';
	time= [time(2), time(2)+tchunk];
	totaltime= totaltime+ tchunk;
	cv= max(var(x)./mean(x));
	%allx= [allx;x];
	%allt= [allt;t];
  if (keep_proportion > 0)
% $$$     keep_indices = fliplr(...
% $$$       floor(length(t):(-1. / keep_proportion):1));
    keep_indices = rand(size(x, 1), 1) <= keep_proportion;
    allx= [allx;x(keep_indices, :)];
    allt= [allt;t(keep_indices, :)];
  end
  stop_condition = ~(cv > maxcv & totaltime < tspan(length(tspan)));
  if (stop_condition && keep_proportion == 0)
    allx = x(end, :); 
    allt = t(end, :); 
  end
	clear t
  %toc
  fprintf('%f time elapsed in %s\n', toc, mfilename);
end
%toc