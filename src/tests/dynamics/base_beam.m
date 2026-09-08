%% Floating rigid base with a rigidly attached flexible beam
% A sinusoidal point force is applied in the global y direction at the tip.
clear; close all; clc;

thisFile = mfilename('fullpath');
testDir = fileparts(thisFile);
srcDir = fileparts(fileparts(testDir));
projectDir = fileparts(srcDir);
dataDir = fullfile(projectDir, 'data', 'base_beam');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

addpath(fullfile(srcDir, 'classes'), fullfile(srcDir, 'utils'));

%% Flexible beam: same properties as flexible_base_manipulator.m
beam.length = 1.0;
beam.width = 0.032;
beam.thickness = 0.0015;
beam.density = 7.542e3;
beam.youngs_modulus = 1.983e11;
beam.poisson_ratio = 0.30;
beam.elements_num = 3;
beam.fixed_root = true;
beam.zeta = 0.0;
beam.omega_n = 7.808;

femBeam = FEM_Beam(beam);
[~, K3d] = femBeam.get_matrices_3d();
freeDofs = femBeam.get_free_dofs();

% Retain [ux, uy, rz] at each free node for planar axial/bending motion.
beam.dof = 3 * beam.elements_num;
beamGlobalDofs = zeros(1, beam.dof);
for node = 2:(beam.elements_num + 1)
    nodeDofs = femBeam.get_node_dofs(node);
    ids = 3 * (node - 2) + (1:3);
    beamGlobalDofs(ids) = nodeDofs([1, 2, 6]);
end
[isFree, beamCoordIds] = ismember(beamGlobalDofs, freeDofs);
if ~all(isFree)
    error('Could not map the retained beam DOFs into the fixed-root matrices.');
end

shapeMatrix = zeros(femBeam.get_total_dof(), beam.dof);
for k = 1:beam.dof
    shapeMatrix(beamGlobalDofs(k), k) = 1;
end

Kff = full(K3d(beamCoordIds, beamCoordIds));
Kff = 0.5 * (Kff + Kff.');
Dff = 2 * beam.zeta / beam.omega_n * Kff;
[u0Beam, Sbeam, beamWeights] = femBeam.build_ffrf_data( ...
    shapeMatrix, 4, 3, 3);

tipShape = femBeam.evaluate_displacement( ...
    shapeMatrix, beam.length, 0, 0);

%% Rigid base included in the same floating frame
% Two-point Gauss integration exactly reproduces the mass and inertia of
% the rectangular base while its elastic shape values remain zero.
base.mass = 5.0;
base.length = 0.20;
base.height = 0.12;
base.depth = 0.08;
base.com = [-0.5 * base.length; 0; 0];
base.volume = base.length * base.height * base.depth;
base.density = base.mass / base.volume;

gaussCoordinate = [-1, 1] / sqrt(3);
[xi, yi, zi] = ndgrid(gaussCoordinate, gaussCoordinate, gaussCoordinate);
u0Base = [base.com(1) + 0.5 * base.length * xi(:).'; ...
          0.5 * base.height * yi(:).'; ...
          0.5 * base.depth * zi(:).'];
nBasePoints = size(u0Base, 2);
Sbase = zeros(3, beam.dof, nBasePoints);
baseWeights = base.volume / nBasePoints * ones(nBasePoints, 1);

compositeParams.u0 = [u0Beam, u0Base];
compositeParams.S = cat(3, Sbeam, Sbase);
compositeParams.weights = [beamWeights(:); baseWeights];
compositeParams.density = [beam.density * ones(numel(beamWeights), 1); ...
                           base.density * ones(nBasePoints, 1)];
compositeParams.Kff = Kff;
compositeParams.Dff = Dff;
compositeParams.q = [1; 0; 0; 0; zeros(3 + beam.dof, 1)];
compositeParams.qdot = zeros(7 + beam.dof, 1);
baseBeamBody = FlexBody(compositeParams);

expectedMass = base.mass + beam.density * beam.length * ...
    beam.width * beam.thickness;
if abs(baseBeamBody.mass - expectedMass) > 1e-10 * expectedMass
    error('Composite base-beam mass does not match the specified geometry.');
end

%% Sinusoidal beam-tip excitation and integration
forceAmplitude = 0.1;   % N
forceFrequency = 1.0;   % Hz
simulationTime = 10.0;   % s
outputStep = 1.0e-3;    % s
tOutput = (0:outputStep:simulationTime).';

nCoordinate = 7 + beam.dof;
x0 = [compositeParams.q; compositeParams.qdot];
odeOptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, ...
    'MaxStep', outputStep);

tic;
[tVec, xHist] = ode15s( ...
    @(t, x) base_beam_rhs(t, x, baseBeamBody, tipShape, ...
        beam.length, forceAmplitude, forceFrequency), ...
    tOutput, x0, odeOptions);
elapsed = toc;

qHist = xHist(:, 1:nCoordinate);
quaternionHist = qHist(:, 1:4);
quaternionHist = quaternionHist ./ vecnorm(quaternionHist, 2, 2);
basePosition = qHist(:, 5:7);
etaHist = qHist(:, 7 + (1:beam.dof));
beamTipElastic = etaHist * tipShape.';

nTime = numel(tVec);
baseTheta = zeros(nTime, 1);
beamTipGlobal = zeros(nTime, 3);
for k = 1:nTime
    A = get_rotation_matrix_quaternion(quaternionHist(k, :).');
    baseTheta(k) = atan2(A(2, 1), A(1, 1));
    localTip = [beam.length; 0; 0] + beamTipElastic(k, :).';
    beamTipGlobal(k, :) = (basePosition(k, :).' + A * localTip).';
end
forceY = forceAmplitude * sin(2 * pi * forceFrequency * tVec);

outputCsv = fullfile(dataDir, 'matlab_base_beam.csv');
writematrix_with_header(outputCsv, ...
    [tVec, baseTheta, beamTipGlobal(:, 2), beamTipElastic(:, 2), ...
     basePosition(:, 1:2), beamTipGlobal(:, 1), forceY], ...
    ['time,base_theta,beam_tip_global_y,beam_tip_elastic_y,', ...
     'base_x,base_y,beam_tip_global_x,force_y']);

%% MATLAB, Exudyn, and Chrono comparison
exudynCsv = fullfile(dataDir, 'exudyn_base_beam.csv');
chronoCsv = fullfile(dataDir, 'chrono_base_beam.csv');
if isfile(exudynCsv) && isfile(chronoCsv)
    exudynData = read_reference_data(exudynCsv, tVec);
    chronoData = read_reference_data(chronoCsv, tVec);

    comparisonFig = figure('Color', 'w', 'Position', [100, 100, 980, 720]);
    layout = tiledlayout(comparisonFig, 2, 1, ...
        'Padding', 'compact', 'TileSpacing', 'compact');
    colors = [0.00, 0.45, 0.74; ...
              0.85, 0.33, 0.10; ...
              0.47, 0.67, 0.19];

    ax = nexttile(layout);
    plot_comparison(ax, tVec, baseTheta, exudynData(:, 1), ...
        chronoData(:, 1), 'Base angle (rad)', colors, true);

    ax = nexttile(layout);
    plot_comparison(ax, tVec, beamTipGlobal(:, 2), exudynData(:, 2), ...
        chronoData(:, 2), 'Beam-tip global y displacement (m)', ...
        colors, false);

    comparisonPng = fullfile(dataDir, 'base_beam_comparison.png');
    exportgraphics(comparisonFig, comparisonPng, 'Resolution', 220);

    exudynRms = sqrt(mean([baseTheta - exudynData(:, 1), ...
        beamTipGlobal(:, 2) - exudynData(:, 2)].^2, 1));
    chronoRms = sqrt(mean([baseTheta - chronoData(:, 1), ...
        beamTipGlobal(:, 2) - chronoData(:, 2)].^2, 1));
    fprintf('Exudyn RMS [base angle, tip y]: [%.4e rad, %.4e m]\n', ...
        exudynRms(1), exudynRms(2));
    fprintf('Chrono RMS [base angle, tip y]: [%.4e rad, %.4e m]\n', ...
        chronoRms(1), chronoRms(2));
    fprintf('Comparison figure: %s\n', comparisonPng);
else
    warning(['Reference data are missing. Run refs/exudyn_base_beam.py ', ...
        'and refs/chrono_base_beam.py, then rerun this script.']);
end

fprintf('MATLAB base-beam simulation completed in %.3f s.\n', elapsed);
fprintf('Data saved to: %s\n', outputCsv);

function dx = base_beam_rhs(t, x, body, tipShape, beamLength, ...
        forceAmplitude, forceFrequency)
nCoordinate = numel(x) / 2;
q = x(1:nCoordinate);
qdot = x(nCoordinate + 1:end);
[q, qdot] = normalize_quaternion_state(q, qdot);

M = body.get_mass_matrix(q);
Qv = body.get_quadratic_velocity_vec(q, qdot);
Qe = body.get_elastic_force(q, qdot);

theta = q(1:4);
eta = q(8:end);
A = get_rotation_matrix_quaternion(theta);
[~, Gbar] = get_quaternion_G_matrices(theta);

forceGlobal = [0; forceAmplitude * sin(2 * pi * forceFrequency * t); 0];
forceLocal = A.' * forceGlobal;
tipLocal = [beamLength; 0; 0] + tipShape * eta;

% Virtual work of the point force at r = R + A*(u_tip + S_tip*eta).
Qext = zeros(nCoordinate, 1);
Qext(1:4) = Gbar.' * cross(tipLocal, forceLocal);
Qext(5:7) = forceGlobal;
Qext(8:end) = tipShape.' * forceLocal;

constraintJacobian = [2 * theta.', zeros(1, nCoordinate - 4)];
constraintAcceleration = -2 * (qdot(1:4).' * qdot(1:4));
solution = [M, -constraintJacobian.'; constraintJacobian, 0] \ ...
    [Qext + Qv - Qe; constraintAcceleration];
dx = [qdot; solution(1:nCoordinate)];
end

function referenceInterp = read_reference_data(csvFile, tVec)
data = readmatrix(csvFile, 'NumHeaderLines', 1);
data = data(all(isfinite(data), 2), :);
data = sortrows(data, 1);
[time, uniqueIds] = unique(data(:, 1), 'stable');
data = data(uniqueIds, :);
referenceInterp = interp1(time, data(:, 2:3), tVec, 'linear', 'extrap');
end

function plot_comparison(ax, t, matlabValue, exudynValue, chronoValue, ...
        yLabelText, colors, showLegend)
plot(ax, t, matlabValue, '-', 'Color', colors(1, :), ...
    'LineWidth', 1.4); hold(ax, 'on');
plot(ax, t, exudynValue, '--', 'Color', colors(2, :), ...
    'LineWidth', 1.3);
plot(ax, t, chronoValue, ':', 'Color', colors(3, :), ...
    'LineWidth', 1.6); hold(ax, 'off');
grid(ax, 'on'); box(ax, 'on');
xlabel(ax, 'Time (s)'); ylabel(ax, yLabelText);
if showLegend
    legend(ax, {'MATLAB FFRF', 'Exudyn', 'Chrono'}, 'Location', 'best');
end
end
