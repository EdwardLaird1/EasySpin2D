function [BList, nu, map, info] = chili2D(Sys, Exp, BList, Opt, doPlot)
%CHILI2D  Simulate a 2D frequency-swept cw EPR map I(B,nu) using EasySpin chili.
%
%   [BList, nu, map] = chili2D(Sys, Exp, BList)
%   [BList, nu, map] = chili2D(Sys, Exp, BList, Opt)
%   [BList, nu, map] = chili2D(Sys, Exp, BList, Opt, doPlot)
%   [BList, nu, map, info] = chili2D(...)
%
% Outputs
%   info    Structure with metadata, including
%             .elapsedTime    elapsed wall-clock time in seconds
%             .nStates        estimated Hilbert-space dimension

  if nargin < 4 || isempty(Opt)
    Opt = struct();
  end
  if nargin < 5 || isempty(doPlot)
    doPlot = true;
  end

  % ---------- Basic validation ----------
  validateattributes(BList, {'numeric'}, ...
    {'vector','nonempty','real','finite'}, mfilename, 'BList', 3);
  BList = BList(:).';  % row vector

  if ~isstruct(Sys)
    error('Sys must be a structure.');
  end
  if ~isstruct(Exp)
    error('Exp must be a structure.');
  end

  hasMwRange = isfield(Exp,'mwRange') && ~isempty(Exp.mwRange);
  hasMwCenterSweep = isfield(Exp,'mwCenterSweep') && ~isempty(Exp.mwCenterSweep);
  if ~(hasMwRange || hasMwCenterSweep)
    error('Exp must contain either Exp.mwRange or Exp.mwCenterSweep for frequency sweeps.');
  end

  if isfield(Exp,'mwFreq') && ~isempty(Exp.mwFreq)
    error('For frequency sweeps, remove Exp.mwFreq and use Exp.Field.');
  end

  if ~isfield(Exp,'nPoints') || isempty(Exp.nPoints)
    Exp.nPoints = 1024;
  end
  if ~isfield(Exp,'Harmonic') || isempty(Exp.Harmonic)
    Exp.Harmonic = 0;
  end

  if ~(isscalar(doPlot) && (islogical(doPlot) || isnumeric(doPlot)))
    error('doPlot must be a scalar logical/numeric value.');
  end
  doPlot = logical(doPlot);

  % ---------- chili-specific validation ----------
  hasDynamics = ...
    (isfield(Sys,'tcorr')    && ~isempty(Sys.tcorr))    || ...
    (isfield(Sys,'logtcorr') && ~isempty(Sys.logtcorr)) || ...
    (isfield(Sys,'Diff')     && ~isempty(Sys.Diff))     || ...
    (isfield(Sys,'logDiff')  && ~isempty(Sys.logDiff));

  if ~hasDynamics
    error(['chili requires a motional parameter in Sys: ', ...
           'tcorr, logtcorr, Diff, or logDiff.']);
  end

  badStrainFields = {'gStrain','AStrain','DStrain','HStrain','lwStrain'};
  for k = 1:numel(badStrainFields)
    if isfield(Sys,badStrainFields{k}) && ~isempty(Sys.(badStrainFields{k}))
      error('chili does not support strains. Remove Sys.%s.', badStrainFields{k});
    end
  end

  if ~hasAnisotropy(Sys)
    error(['This is an isotropic spin system. chili2D only works for anisotropic ', ...
           'systems, just like chili itself.']);
  end

  % ---------- Timer and metadata ----------
  tStart = tic;
  nB = numel(BList);
  nStates = estimateSpinStates(Sys);

  % ---------- First simulation ----------
  Exp1 = Exp;
  Exp1.Field = BList(1);
  [nu, spc1, infoFirst] = chili(Sys, Exp1, Opt);

  nu = nu(:).';
  spc1 = spc1(:).';

  map = zeros(nB, numel(nu));
  map(1,:) = spc1;

  tol = 1e-10 * max(1, max(abs(nu)));

  % ---------- Remaining fields ----------
  for k = 2:nB
    Expk = Exp;
    Expk.Field = BList(k);

    [nuk, spck] = chili(Sys, Expk, Opt);

    nuk = nuk(:).';
    spck = spck(:).';

    if numel(nuk) ~= numel(nu) || any(abs(nuk - nu) > tol)
      error(['Frequency axis changed between simulations at B = %.6g mT. ', ...
             'Use identical Exp.mwRange/Exp.mwCenterSweep and Exp.nPoints.'], BList(k));
    end

    map(k,:) = spck;
  end

  elapsedTime = toc(tStart);

  % ---------- Optional info output ----------
  if nargout >= 4
    info = struct();
    info.infoFirst = infoFirst;
    info.nFields = nB;
    info.nFreq = numel(nu);
    info.BList = BList;
    info.nu = nu;
    info.elapsedTime = elapsedTime;
    info.nStates = nStates;
  end

  % ---------- Optional plot ----------
  if doPlot
    scale = 0.02 * max(abs(map(:)));
    if ~isfinite(scale) || scale == 0
      scale = 1;
    end

    imagesc(BList, nu, asinh(map/scale).')
    axis xy tight
    xlabel('Magnetic field (mT)')
    ylabel('Microwave frequency (GHz)')

    if isfield(Exp,'mwMode') && ~isempty(Exp.mwMode)
      modeLabel = char(lower(string(Exp.mwMode)));
      modeLabel(1) = upper(modeLabel(1));
      title(sprintf('%s-mode 2D chili cw spectral map', modeLabel))
    else
      title('2D chili cw spectral map')
    end

    cb = colorbar;
    cb.Label.String = 'asinh(signal / scale)';
    box on

    if nargout >= 4
      info.plotScale = scale;
    end
  end

  fprintf('Chili2D: completed in %.3f s for %d x %d points in an estimated %d-state system.\n', ...
    elapsedTime, nB, numel(nu), nStates);
end

function tf = hasAnisotropy(Sys)
% Return true if common spin-Hamiltonian terms are anisotropic.

  tf = false;
  tol = 1e-12;

  if isfield(Sys,'g') && ~isempty(Sys.g)
    g = Sys.g;
    if numel(g) > 1 && any(abs(g(:) - g(1)) > tol)
      tf = true;
      return
    end
  end

  if isfield(Sys,'A') && ~isempty(Sys.A)
    A = Sys.A;
    if size(A,2) >= 3
      for k = 1:size(A,1)
        row = A(k,:);
        if any(abs(row - row(1)) > tol)
          tf = true;
          return
        end
      end
    end
  end

  if isfield(Sys,'D') && ~isempty(Sys.D)
    D = Sys.D;
    if numel(D) > 1 && any(abs(D(:) - D(1)) > tol)
      tf = true;
      return
    end
  end

  if isfield(Sys,'Q') && ~isempty(Sys.Q)
    Q = Sys.Q;
    if size(Q,2) >= 3
      for k = 1:size(Q,1)
        row = Q(k,:);
        if any(abs(row - row(1)) > tol)
          tf = true;
          return
        end
      end
    end
  end
end

function nStates = estimateSpinStates(Sys)
% Estimate Hilbert-space dimension from electron and nuclear spins.

  if isfield(Sys,'S') && ~isempty(Sys.S)
    if numel(Sys.S) > 1
      nStates = prod(2*Sys.S(:) + 1);
    else
      nStates = 2*Sys.S + 1;
    end
  else
    nStates = 1;
  end

  if isfield(Sys,'Nucs') && ~isempty(Sys.Nucs)
    nucs = strtrim(split(string(Sys.Nucs), ','));
    nucs = nucs(nucs ~= "");

    if isfield(Sys,'n') && ~isempty(Sys.n)
      mults = Sys.n(:).';
      if numel(mults) ~= numel(nucs)
        error('Sys.n must have the same number of elements as nuclei in Sys.Nucs.');
      end
    else
      mults = ones(1,numel(nucs));
    end

    for k = 1:numel(nucs)
      I = nucspin(char(nucs(k)));
      nStates = nStates * (2*I + 1)^mults(k);
    end
  end

  nStates = round(nStates);
end