%% Demo: garlic2D, pepper2D, and chili2D
% Open this script in MATLAB Live Editor if you want a live-script view,
% then save as .mlx.
%
% Place this file in the same folder as garlic2D.m, pepper2D.m, and chili2D.m,
% or add those folders to the MATLAB path.

clearvars
close all
clc

thisFile = mfilename('fullpath');
if ~isempty(thisFile)
    addpath(fileparts(thisFile));
end

%% Check EasySpin
% This prints the installed EasySpin version.
easyspin

%% Define a minimal anisotropic S = 1/2 test system
% This system is simple enough for a compact demo but anisotropic enough
% to exercise chili2D and pepper2D.

Radical = 'Minimal anisotropic S = 1/2 demo';

Sys = struct();
Sys.S = 1/2;
Sys.g = [2.010 2.006 2.002];      % principal values
Sys.Nucs = '14N';
Sys.A = [20 20 95];               % MHz

isAnisotropic = hasAnisotropy(Sys);

%% Experimental settings for frequency-swept 2D maps
Exp = struct();
Exp.mwRange = [8.5 10.5];         % GHz
Exp.nPoints = 400;
Exp.Harmonic = 0;

BList = linspace(300,360,121);    % mT

Opt = struct();
Opt.Verbosity = 0;

%% Garlic2D: fast-motion solution regime
% Garlic is intended for S = 1/2 radicals in the isotropic / fast-motion
% regime. Here we add a short correlation time.

SysGarlic = Sys;
SysGarlic.tcorr = 3e-11;          % s, fast motion
SysGarlic.lwpp = 0.15;            % linewidth for garlic demo

Exp.mwMode = 'perpendicular';
[BgPerp, nuGPerp, mapGPerp, infoGPerp] = garlic2D(SysGarlic,Exp,BList,Opt,false);

Exp.mwMode = 'parallel';
[BgPar, nuGPar, mapGPar, infoGPar] = garlic2D(SysGarlic,Exp,BList,Opt,false);

plotTwoMaps(BgPerp,nuGPerp,mapGPerp, ...
            BgPar,nuGPar,mapGPar, ...
            sprintf('garlic2D — %s', Radical), ...
            'scale_G');

%% Pepper2D: rigid-limit / solid-state regime
SysPepper = Sys;
SysPepper.lwpp = 0.15;

Exp.mwMode = 'perpendicular';
[BpPerp, nuPPerp, mapPPerp, infoPPerp] = pepper2D(SysPepper,Exp,BList,Opt,false);

Exp.mwMode = 'parallel';
[BpPar, nuPPar, mapPPar, infoPPar] = pepper2D(SysPepper,Exp,BList,Opt,false);

plotTwoMaps(BpPerp,nuPPerp,mapPPerp, ...
            BpPar,nuPPar,mapPPar, ...
            sprintf('pepper2D — %s', Radical), ...
            'scale_P');

%% Chili2D: slow-motion regime (perpendicular mode only)
% Current EasySpin chili does not support parallel-mode spectra.

if isAnisotropic
    SysChili = Sys;
    SysChili.tcorr = 3e-9;                    % s, slow tumbling
    SysChili.lwpp = unitconvert(0.15,'mT->MHz');

    Exp.mwMode = 'perpendicular';
    [BcPerp, nuCPerp, mapCPerp, infoCPerp] = chili2D(SysChili,Exp,BList,Opt,false);

    plotOneMap(BcPerp,nuCPerp,mapCPerp, ...
               sprintf('chili2D — %s', Radical), ...
               'Perpendicular mode only', ...
               'scale_C');
else
    disp('Skipping chili2D demo: system is isotropic.')
end

%% Summary
fprintf('\nDemo complete.\n');
fprintf('garlic2D points: %d x %d\n', infoGPerp.nFields, infoGPerp.nFreq);
fprintf('pepper2D points: %d x %d\n', infoPPerp.nFields, infoPPerp.nFreq);
if isAnisotropic
    fprintf('chili2D points:  %d x %d\n', infoCPerp.nFields, infoCPerp.nFreq);
end

%% Local functions
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

function plotTwoMaps(B1,nu1,map1,B2,nu2,map2,figTitle,colorbarLabel)
% Plot two 2D maps side by side with a shared top colorbar.

scale = 0.02*max(abs([map1(:); map2(:)]));
if ~isfinite(scale) || scale==0
    scale = 1;
end

disp1 = asinh(map1/scale);
disp2 = asinh(map2/scale);
clims = [min([disp1(:); disp2(:)]), max([disp1(:); disp2(:)])];

fig = figure('Position',[100 100 1000 700]);
t = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(t);
imagesc(ax1,B1,nu1,disp1.')
axis(ax1,'xy','tight')
xlabel(ax1,'Magnetic field (mT)')
ylabel(ax1,'Microwave frequency (GHz)')
title(ax1,{figTitle,'Perpendicular mode'})
clim(ax1,clims)
box(ax1,'on')

ax2 = nexttile(t);
imagesc(ax2,B2,nu2,disp2.')
axis(ax2,'xy','tight')
xlabel(ax2,'Magnetic field (mT)')
ylabel(ax2,'Microwave frequency (GHz)')
title(ax2,{figTitle,'Parallel mode'})
clim(ax2,clims)
box(ax2,'on')

colormap(fig,parula)
cb = colorbar;
cb.Layout.Tile = 'north';
cb.Label.String = sprintf('asinh(signal / %s)', colorbarLabel);

linkaxes([ax1 ax2],'xy')
end

function plotOneMap(B,nu,map,figTitle,subtitleText,colorbarLabel)
% Plot one 2D map with a top colorbar.

scale = 0.02*max(abs(map(:)));
if ~isfinite(scale) || scale==0
    scale = 1;
end

dispMap = asinh(map/scale);
clims = [min(dispMap(:)), max(dispMap(:))];

fig = figure('Position',[150 150 700 700]);
ax = axes(fig);
imagesc(ax,B,nu,dispMap.')
axis(ax,'xy','tight')
xlabel(ax,'Magnetic field (mT)')
ylabel(ax,'Microwave frequency (GHz)')
title(ax,{figTitle,subtitleText})
clim(ax,clims)
box(ax,'on')

colormap(fig,parula)
cb = colorbar(ax,'Location','northoutside');
cb.Label.String = sprintf('asinh(signal / %s)', colorbarLabel);
end
