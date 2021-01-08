function dfdf = fp_jacobian(t,y,params,adjoint)
if nargin < 4
  adjoint = false;
end
switch params.protocol
case 'Current'
  if isstruct(y)
    y = sol_interp(y,t);
  end
  V = y(end);
  y = y(1:end-1);
case 'Voltage'
  V = params.V(t);
end
dc = params.dc;
n = params.nc;
if isstruct(params.OCV)
  OCV = params.OCV.val;
else
  OCV = params.OCVval;
end
eta = V - OCV;
R = BVSR(params.R0.val,eta,params.alpha,params.RS,params.unit);
%R./(-eta)
Reta = R./(-eta);
ind = (eta==0);
Reta(ind) = params.R0.val(ind) ./ (1+params.R0.val(ind)*params.RS) ./params.const;
D = (params.D0).*Reta;

R_extend=[0;R;0];
dfdf_lower=[R;0];
dfdf_upper=-[0;R];
dfdf=[dfdf_lower,-diff(R_extend,1,1),dfdf_upper]/2/dc ...
    +[[D;0],[-D(1);-D(1:end-1)-D(2:end);-D(end)],[0;D]]/dc^2;
dfdf=spdiags(dfdf,[-1,0,1],n,n);

if isequal(params.protocol,'Current')
  %differentiation of fokker_planck source term with respect to voltage
  dR = BVSRdeta(params.R0.val,eta,params.alpha,params.RS,R,params.unit);
  dD = (params.D0).*(dR./(-eta) + R./eta.^2);
  F = dR .* (y(1:end-1)+y(2:end))/2;
  F = F - dD.*diff(y)/dc;
  F = [0;F;0];
  dy = -diff(F)/dc;
  dfdf = [dfdf,dy];
  dfdf(n+1,n+1) = 0;
end

if adjoint
  dfdf = dfdf';
end
