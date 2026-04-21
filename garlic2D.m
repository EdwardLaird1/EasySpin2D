function [BList, nu, map, info] = garlic2D(Sys, Exp, BList, Opt, doPlot)
%GARLIC2D  Simulate a 2D frequency-swept cw EPR map I(B,nu) using EasySpin garlic.
%
%   [BList, nu, map] = garlic2D(Sys, Exp, BList)
%   [BList, nu, map] = garlic2D(Sys, Exp, BList, Opt)
%   [BList, nu, map] = garlic2D(Sys, Exp, BList, Opt, doPlot)
%   [BList, nu, map, info] = garlic2D(...)
%
% Inputs
%   Sys     EasySpin spin system structure for garlic.
%           Intended for the same regime as garlic itself:
%           isotropic / fast-motion solution radicals, typically S = 1/2.
%
%   Exp     EasySpin experiment structure for FREQUENCY-SWEPT simulations.
%           Required:
%             Exp.mwRange   or Exp.mwCenterSweep
%           Must NOT contain:
%             Exp.mwFreq
%           Should NOT contain:
%             Exp.ModAmp
%
%   BList   Vector of static magnetic fields in mT
%
%   Opt     EasySpin options structure (optional)
%
%   doPlot  true/false, whether to plot the result (default = true)
%
% Outputs
%   BList   Field axis in mT (row vector)
%   nu      Frequency axis in GHz (row vector)
%   map     2D matrix of size [numel(BList) x numel(nu)]
%           containing the simulated frequency-swept cw spectral intensity
%   info    Structure with metadata:
%             .infoFirst      info output from the first garlic call
%             .nFields        number of field points
%             .nFreq          number of frequency points
%             .BList          field axis (mT)
%             .nu             frequency axis (GHz)
%             .plotScale      scale used in asinh plotting
%             .elapsedTime    elapsed wall-clock time in seconds
%             .nStates        estimated Hilbert-space dimension
%
% Notes
%   - This function returns a stacked-spectra 2D cw map, not a discrete
%     transition-trajectory map.
%   - Exp.ModAmp is not supported here for frequency sweeps.
%   - For fast-motion simulations (Sys.tcorr or Sys.logtcorr nonzero),
%     Opt.AccumMethod = 'binning' and 'template' are rejected because
%     garlic does not support them in that regime.

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

  if ~isstruct(Exp)
    error('Exp must be a structure.');
  end

  hasMwRange = isfield(Exp,'mwRange');
  hasMwCenterSweep = isfield(Exp,'mwCenterSweep');
  if ~(hasMwRange || hasMwCenterSweep)
    error('Exp must contain either Exp.mwRange or Exp.mwCenterSweep for frequency sweeps.');
  end

  if isfield(Exp,'mwFreq')
    error('For frequency sweeps, remove Exp.mwFreq and use Exp.Field.');
  end

  if isfield(Exp,'ModAmp') && any(Exp.ModAmp ~= 0)
    error('Exp.ModAmp is not supported in this wrapper for frequency sweeps.');
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

  % ---------- Lightweight "garlic-like" sanity check ----------
  % garlic is intended for S = 1/2 solution radicals.
  if ~iscell(Sys) && isfield(Sys,'S') && ~isempty(Sys.S)
    if ~(isscalar(Sys.S) && abs(Sys.S - 1/2) < 1e-12)
      warning(['Sys.S = %g. garlic is intended for solution radicals with S = 1/2; ', ...
               'results may be invalid for this system.'], Sys.S);
    end
  end

  % ---------- Fast-motion guard for unsupported accumulation methods ----------
  hasFastMotion = false;
  if ~iscell(Sys)
    if isfield(Sys,'logtcorr') && ~isempty(Sys.logtcorr)
      hasFastMotion = true;
    elseif isfield(Sys,'tcorr') && ~isempty(Sys.tcorr) && any(Sys.tcorr ~= 0)
      hasFastMotion = true;
    end
  end

  if hasFastMotion && isfield(Opt,'AccumMethod') && ~isempty(Opt.AccumMethod)
    badAccum = any(strcmpi(Opt.AccumMethod, {'binning','template'}));
    if badAccum
      error(['Opt.AccumMethod=''%s'' is not valid for fast-motion garlic simulations. ', ...
             'Use ''explicit'' or leave it empty.'], Opt.AccumMethod);
    end
  end

  % ---------- Timer and metadata ----------
  tStart = tic;
  nB = numel(BList);
  nStates = estimateSpinStates(Sys);

  % ---------- First simulation: establish common frequency axis ----------
  Exp1 = Exp;
  Exp1.Field = BList(1);
  [nu, spc1, infoFirst] = garlic(Sys, Exp1, Opt);

  nu = nu(:).';
  spc1 = spc1(:).';

  map = zeros(nB, numel(nu));
  map(1,:) = spc1;

  % Relative numerical tolerance for axis consistency
  tol = 1e-10 * max(1, max(abs(nu)));

  % ---------- Remaining fields ----------
  for k = 2:nB
    Expk = Exp;
    Expk.Field = BList(k);

    [nuk, spck] = garlic(Sys, Expk, Opt);

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
      title(sprintf('%s-mode 2D garlic cw spectral map', modeLabel))
    else
      title('2D garlic cw spectral map')
    end

    cb = colorbar;
    cb.Label.String = 'asinh(signal / scale)';
    box on

    if nargout >= 4
      info.plotScale = scale;
    end
  end

  fprintf('Garlic2D: completed in %.1f s for %d x %d points in an estimated %d-state system.\n', ...
    elapsedTime, nB, numel(nu), nStates);
end

function nStates = estimateSpinStates(Sys)
% Estimate Hilbert-space dimension from electron and nuclear spins.

  % Electron-spin contribution
  if isfield(Sys,'S') && ~isempty(Sys.S)
    if numel(Sys.S) > 1
      nStates = prod(2*Sys.S(:) + 1);
    else
      nStates = 2*Sys.S + 1;
    end
  else
    nStates = 1;
  end

  % Nuclear-spin contribution
  if isfield(Sys,'Nucs') && ~isempty(Sys.Nucs)
    nucs = strtrim(split(string(Sys.Nucs), ','));
    nucs = nucs(nucs ~= "");

    % Multiplicities for equivalent nuclei
    if isfield(Sys,'n') && ~isempty(Sys.n)
      mults = Sys.n(:).';
      if numel(mults) ~= numel(nucs)
        error('Sys.n must have the same number of elements as nuclei in Sys.Nucs.');
      end
    else
      mults = ones(1,numel(nucs));
    end

    for k = 1:numel(nucs)
      I = nucspin(char(nucs(k)));   % EasySpin nuclear spin
      nStates = nStates * (2*I + 1)^mults(k);
    end
  end

  nStates = round(nStates);
end