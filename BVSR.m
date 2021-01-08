function I = BVSR(I0,eta,alpha,RS,unit,AbsTol,RelTol)
%Butler-Volmer with series resistance, the current is given by
%I0*(exp(alpha*(-eta-RS*I))-exp(-(1-alpha)*(-eta-RS*I))) = I
%unit specifies the unit of eta, by default it is 'kbT', it can also be 'V' (volt)
%RS is defined in kbT/I0 unit (whatever unit I0 has)
%I0 and eta can be arrays of compatible sizes
%for example, I can be computed for an array of c at the same voltage, R0 should then contain exchange current at the requested c and eta should be the overpotential at those c values.
%I can also be computed for an array of voltage at the same c, in which case R0 should be a scalar and eta should correspond to the overpotential at the same c.
%Or eta can be constant, in which case we compute the current/reaction at different c values with the same overpotential.
%Or R0 can be an nc*1 array and eta can be a nc*nt matrix.
%AbsTol and RelTol are absolute and relative tolerance on the convergence criteria based on the change in I after each Newton iteration
if nargin < 3 || isempty(alpha)
  alpha = 0.5;
end
if nargin < 4 || isempty(RS)
  RS = 0;
end
if nargin < 5 || isempty(unit)
  unit = 'kbT';
end
if nargin < 6 || isempty(AbsTol)
  AbsTol = 1e-8;
end
if nargin < 7 || isemty(RelTol)
  RelTol = 1e-6;
end

switch unit
case 'kbT'
  const = 1;
case 'V'
  const = 0.025693; %RT/F 25.693mV
end
eta = eta/const;

Ig1 = I0.*(exp(alpha.*(-eta))-exp(-(1-alpha).*(-eta)));
if RS == 0
  I = Ig1;
  return;
else
  %another guess, in the limit of RS being very large, I \approx -eta/RS
  Ig2 = -eta./RS;
  I = zeros(size(Ig1));
  ind = abs(Ig1)<0.01*abs(Ig2);
  I(ind) = Ig1(ind);
end

Isize = size(I);
% Newton's method
while true
   Iprev = I;
   lnp = alpha.*(-eta-I.*RS);
   lnm = -(1-alpha).*(-eta-I.*RS);
   scale = max(lnp,lnm);
   p = exp(lnp-scale);
   m = exp(lnm-scale);
   f = I0.*(p-m)-I.*exp(-scale);
   df = I0.*(-alpha.*RS.*p-(1-alpha).*RS.*m)-exp(-scale);
   I = I-f./df;
   if any(isnan(I(:)))
     %once NaN shows up, it will never go away so it'll be stuck in
     %the while loop. Abort and set all entries of I to NaN
     I = NaN(Isize);
     break;
   elseif all(abs(I - Iprev) < AbsTol + RelTol*abs(Iprev))
     break;
   end
end
