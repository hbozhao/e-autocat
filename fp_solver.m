function [tout,y,params,mass] = fp_solver(x,tspan,y0,params,solver,sol)
%Fokker-Planck equation solver
%x is the mesh (uniform), tspan is requested time point, y0 is the initial condition defined on x
%supply the following parameters and functions in params struct
%params.D0 is the D0 in the diffusive term
%params.OCV can either be a function handle for OCV, the input is concentration x
%or a struct which contains func, the function handle and params, the parameters. The usage is func(x,params)
%params.R0 is a struct for exchange current. It may contain function handles, func for function and sensitivity evaluations and params. Alternatively, it may contain precomputed values stored in val (must be evaluated at x_mid), otherwise they will be evaluated inside this function with the function handles and params.
%The former is recommended if exponentiation of the params is needed (as is often the case in inverse problems)
%the usage of func is func(x,params)
%params.V is the dependence of voltage on time. The usage is similar to OCV,
%it can either be a function handle with input time, or a struct that contains func and params
%params.R is a function handle for current, the input is time
%params.protocol can be Current or Voltage, which determines whether to control the voltage using params.V or current using params.R
%optional:
%params.unit is the unit of voltage (kbT, or V), kbT by default
%params.options can be used to set ODE options before fp_solver makes mandatory changes
%use solver to choose what ODE solver in matlab to use. Force to be ode15s if params.protocol = Current
%set sol = true to output sol struct (see ode solvers). false by default (in which which the output y is the solution array)
%the output tout is time returned by the ODE solver, y is either the solution array or sol struct returned by the odesolver, mass is the mass matrix of the ODE
if nargin<5 || isempty(solver)
  solver = 'ode15s';
end
if nargin<6 || isempty(sol)
  sol = false;
end
if size(x,1)==1
  x = x';
end
if size(y0,1)==1
  y0 = y0';
end
x_mid=(x(1:end-1)+x(2:end))/2;

if ~isfield(params,'unit') || isempty(params.unit)
  params.unit = 'kbT';
end
switch params.unit
case 'kbT'
  params.const = 1;
case 'V'
  params.const = 0.025693; %RT/F 25.693mV
end

if ~isfield(params.R0,'val')
  params.R0.val = params.R0.func(x_mid,params.R0.params);
end
if isstruct(params.OCV)
  %OCV is a struct just like R0
  %must be in this form if sensivitity wrt OCV is required
  params.OCV.val = params.OCV.func(x_mid,params.OCV.params);
else
  params.OCVval = params.OCV(x_mid);
end

params.dc=x(2)-x(1);

n = length(x);
params.nc = n; %number of c variables
options = odeset;
if isfield(params,'options') && isstruct(params.options)
  options = odeset(options,params.options);
end
if isequal(params.protocol,'Current')
  mass = spdiags([ones(n,1);0],0,n+1,n+1);
  mass(end,1:n) = x*params.dc;
  options = odeset(options,'Mass',mass,'MassSingular','yes','MStateDependence','none');
  solver = 'ode15s';
  %consistent initialization
  t0 = tspan(1);
  mass_last_row = mass(end,:);
  if length(y0)>n
    dphi0 = y0(n+1);
    y0 = y0(1:n);
  else
    cavg0 = sum(y0(:).*x(:)) / sum(y0);
    if isstruct(params.OCV)
      dphi0 = params.OCV.func(cavg0,params.OCV.params);
    else
      dphi0 = params.OCV(cavg0);
    end
  end
  if ~isfinite(dphi0)
    tout = t0;
    y = [y0(:).',dphi0];
    mass = [];
    return;
  end
  if isstruct(params.R)
    current0 = params.R.func(t0,params.R.params);
  else
    current0 = params.R(t0);
  end
  dphi0 = fzero(@(dphi) mass_last_row*fokker_planck(t0,[y0;dphi],params)-current0,dphi0);
  y0(n+1) = dphi0;
  yp0 = fokker_planck(t0,y0,params);
  yp0(end) = 0;
  options = odeset(options,'InitialSlope',yp0);
else
  mass = speye(n);
end

if isequal(solver,'ode15s')
  options = odeset(options,'Jacobian',@(t,y) fp_jacobian(t,y,params));
end

odeFcn = @(t,y) fokker_planck(t,y,params);
if sol
  y = feval(solver,odeFcn,tspan,y0,options);
  tout = y.x;
else
  [tout,y] = feval(solver,odeFcn,tspan,y0,options);
end
end

function dy=fokker_planck(t,y,params)
dc = params.dc;
switch params.protocol
case 'Current'
  V = y(end);
  y = y(1:end-1);
case 'Voltage'
  V = params.V(t);
end
if isstruct(params.OCV)
  OCV = params.OCV.val;
else
  OCV = params.OCVval;
end
eta = V - OCV;
R = BVSR(params.R0.val,eta,params.alpha,params.RS,params.unit);
if isequal(params.D0mode,'pd')
  %R./(-eta)
  Reta = R./(-eta);
  ind = (eta==0);
  Reta(ind) = params.R0.val(ind) ./ (1+params.R0.val(ind)*params.RS) ./params.const;
end
D = (params.D0).*Reta;
F = R .* (y(1:end-1)+y(2:end))/2;
F = F - D.*diff(y)/dc;
F = [0;F;0];
dy = -diff(F)/dc;
if isequal(params.protocol,'Current')
  if isstruct(params.R)
    dy(end+1) = params.R.func(t,params.R.params);
  else
    dy(end+1) = params.R(t);
  end
end
end
