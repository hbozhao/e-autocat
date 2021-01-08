function yeq = equilibrate_avg(OCV,x,cavg,D0)
  %computes the equilibrium distribution at the desired average concentration cavg.
  %This can be used directly as the initial condition for pd
  %OCV usage is explained in fp_solver
  %x is the mesh to compute the distribution at
  %the output yeq concatenates the distribution and the voltage V i
  %at equilibrium, the distribution is proportional to
  %exp(-\int{(V-OCV)dc}/D0)
  if isstruct(OCV)
    %OCV is a struct just like R0
    %must be in this form if sensivitity wrt OCV is required
    OCVval = OCV.func(x,OCV.params);
    Vguess = OCV.func(cavg,OCV.params);
  else
    OCVval = OCV(x);
    Vguess = OCV(cavg);
  end

  V = fzero(@(V) avg(x,eqpdf(OCvalV,x,V,D0))-cavg,Vguess);
  yeq = eqpdf(OCV,x,V,D0);
  yeq(end+1) = V;
end

function y = eqpdf(OCV,x,V,D0)
  G = cumtrapz(x,V-OCV);
  G = G-min(G);
  y = exp(-G/D0);
  y = y / trapz(x, y);
end

function y=avg(mid,f)
%calculates the mean of X, i.e. \int{xf(x)dx}
%mid is the value of X in the middle of each cell and the f is the corresponding probability density
%assume uniform grid
y = trapz(x, x.*y);
end
