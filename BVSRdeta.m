function dIdeta = BVSRdeta(I0,eta,alpha,RS,I,unit,varargin)
  %The derivative of the current I from the model below (Butler-Volmer with series resistance) with respect to eta
  %I0*(exp(alpha*(-eta-RS*I))-exp(-(1-alpha)*(-eta-RS*I))) = I
  %I is the solution to BVSR, which can be provided to accelerate the evaluation here
  %varargin are for additional arguments of BVSR
  if nargin < 3 || isempty(alpha)
    alpha = 0.5;
  end
  if nargin < 4 || isempty(RS)
    RS = 0;
  end
  if nargin < 6 || isempty(unit)
    unit = 'kbT';
  end
  switch unit
  case 'kbT'
    const = 1;
  case 'V'
    const = 0.025693; %RT/F 25.693mV
  end
  eta = eta/const;

  if RS == 0
    dIdeta = I0.*(-alpha.*exp(alpha.*(-eta))-(1-alpha).*exp(-(1-alpha).*(-eta)));
  else
    if nargin<5
      I = BVSR(I0,eta,alpha,RS,'kbT',varargin{:});
    end
    lnp = alpha.*(-eta-I.*RS);
    lnm = -(1-alpha).*(-eta-I.*RS);
    scale = max(lnp,lnm);
    p = exp(lnp-scale);
    m = exp(lnm-scale);
    %d(equation)/dI
    dfdI = I0.*(-alpha.*RS.*p-(1-alpha).*RS.*m)-exp(-scale);
    %d(equation)/deta
    dfdeta = I0.*(-alpha.*p-(1-alpha).*m);
    dIdeta = -dfdeta./dfdI;
  end
  dIdeta = dIdeta/const;
