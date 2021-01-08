function [PS,cavg,x] = PS_sim(params,I,crange,varargin)
  %computes the phase diagram of instability, that is, current and average concentration at which bimodal distribution occurs
  %the particle ensemble reacts at constant total reaction rate (current)
  %I is the list of current to sweep over
  %crange is the range of average concentration that the ensemble reacts in (a two-element array)
  %see below for varargin options
  %the output PS (phase diagram values) is the prominence value of the second peak, 0 when the distribution is unimodal

  ps = inputParser;
  ps.KeepUnmatched = true;
  addParameter(ps,'OutputNp',100); %the number of point in time (or in terms of average concentration)
  addParameter(ps,'pd_x',linspace(0,1,1000)); %grid point for pd model
  addParameter(ps,'pd_sigma0',[]); %if empty, the initial pdf is fully equilibrated, if a scalar is provided, it is the sigma of the normal pdf used as the initial condition for equilibration for pd
  addParameter(ps,'MinSeparation',0.1); %minimum separation between two peaks in terms of soc
  addParameter(ps,'overwrite_params',[]); %a struct that overwrite existing params,
  %the following will still be forced to change in this function: protocol, R

  parse(ps,varargin{:});
  ps = ps.Results;


  if ~isempty(ps.overwrite_params)
    names = fieldnames(ps.overwrite_params);
    for i = 1:length(names)
      params.(names{i}) = ps.overwrite_params.(names{i});
    end
  end

  OutputNp = ps.OutputNp;
  crange = sort(crange);
  cavg = linspace(crange(1),crange(2),OutputNp);


  params.protocol = 'Current';
  %equilibration
  x = ps.pd_x;
  Nx = length(x);
  f0 = zeros(length(x)+1,2);
  for i = 1:2
    if isempty(ps.pd_sigma0)
      %use fully equilibrated pdf
      f0(:,i) = equilibrate_avg(params.OCV,x,crange(i),params.D0);
    else
      y0 = normpdf(x,crange(i),ps.pd_sigma0);
      cavg0 = sum(y0(:).*x(:))./sum(y0);
      params.R = @(t) t.*exp(-t)*(crange(i)-cavg0); %offset due to error
      tspan = [0,20,50];
      [~,f] = fp_solver(x,tspan,y0,params);
      f0(:,i) = f(end,:);
    end
  end

  for i = 1:length(I)
    thisparams = params;
    thisparams.R = @(t) I(i);
    tf = diff(crange)/abs(I(i));
    tspan = linspace(0,tf,OutputNp);
    [~,y] = fp_solver(x,tspan,f0(:,(I(i)<0)+1),thisparams);
    ff = y(:,1:end-1);
    dc = x(2)-x(1);

    %determine when a second peak forms (defined as phase separation, PS for now)
    %change output to prominence
    %PS is now the second large prominence
    [TF,P] = islocalmax(ff',1,'MaxNumExtrema',2,'MinSeparation',ps.MinSeparation/dc);
    p2 = (sum(TF,1)>1); %columns with 2 extrema
    thisPS = zeros(1,length(p2));
    thisPS(p2) = min(reshape(P(p2&TF),2,[]),[],1);
    if I(i)<0
      thisPS = flip(thisPS);
    end
    PS(i,:) = thisPS;
  end
