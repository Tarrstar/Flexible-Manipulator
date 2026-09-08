%% Floating rigid base connected to a flexible beam by a revolute joint
clear; close all; clc;

thisFile = mfilename('fullpath');
testDir = fileparts(thisFile);
srcDir = fileparts(fileparts(testDir));
projectDir = fileparts(srcDir);
dataDir = fullfile(projectDir, 'data', 'hinged_floating_beam');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

addpath(fullfile(srcDir, 'dynamics'), fullfile(srcDir, 'classes'), ...
    fullfile(srcDir, 'utils'), fullfile(testDir, 'support'));

%% Flexible beam: same parameters as flexible_base_manipulator.m
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

beam.dof = 3 * beam.elements_num;
beamGlobalDofs = zeros(1, beam.dof);
for node = 2:(beam.elements_num + 1)
    nodeDofs = femBeam.get_node_dofs(node);
    ids = 3 * (node - 2) + (1:3);
    beamGlobalDofs(ids) = nodeDofs([1, 2, 6]);
end
[isFree, beamCoordIds] = ismember(beamGlobalDofs, freeDofs);
if ~all(isFree)
    error('Could not map the retained beam DOFs into the fixed-root FE matrices.');
end

shapeMatrix = zeros(femBeam.get_total_dof(), beam.dof);
for k = 1:beam.dof
    shapeMatrix(beamGlobalDofs(k), k) = 1;
end

Kff = full(K3d(beamCoordIds, beamCoordIds));
Dff = 2 * beam.zeta / beam.omega_n * Kff;
[u0, Sbeam, weights] = femBeam.build_ffrf_data(shapeMatrix, 4, 3, 3);

beamParams.density = beam.density;
beamParams.u0 = u0;
beamParams.S = Sbeam;
beamParams.weights = weights;
beamParams.Kff = Kff;
beamParams.Dff = Dff;
beamParams.q = [1; 0; 0; 0; zeros(3 + beam.dof, 1)];
beamParams.qdot = zeros(7 + beam.dof, 1);
flexibleBeam = FlexBody(beamParams);

tipNodeDofs = femBeam.get_node_dofs(beam.elements_num + 1);
tipShape6d = shapeMatrix(tipNodeDofs, :);
beamTipInterface = [tipShape6d(4:6, :); tipShape6d(1:3, :)];
beamTipDisplacement = femBeam.evaluate_displacement( ...
    shapeMatrix, beam.length, 0, 0);

%% Five-kilogram floating rigid base
base.mass = 5.0;
base.length = 0.20;
base.height = 0.12;
base.depth = 0.08;
base.com = [-0.5 * base.length; 0; 0];
base.inertia = base.mass / 12 * [ ...
    base.height^2 + base.depth^2; ...
    base.length^2 + base.depth^2; ...
    base.length^2 + base.height^2];
baseBody = RigidBody(base);

%% Floating base and hinged beam
NB = 2;
model.NB = NB;
model.parent = [0; 1];
model.jtype = {'Rz'; 'Rz'};
% Link 1 is the floating-base placeholder; link 2 is the active hinge.
model.fixed_joints = 1;
model.Xtree = {eye(6); eye(6)};
model.flexbody = [{[]}; {flexibleBeam}];
model.rigidbody = [{baseBody}; {[]}];
model.link_type = {'rigid'; 'flexible'};
model.Phi = [{zeros(6, 0)}; {beamTipInterface}];
model.gravity = zeros(6, 1);
model.floating_base = true;

%% Integration
simulationTime = 2.0;
outputStep = 1.0e-3;
hingeTorque = 0.02;
tOutput = (0:outputStep:simulationTime).';

% x = [base quaternion(4); base position(3); q; eta;
%      base spatial velocity(6); qDot; etaDot].
x0 = zeros(13 + 2 * NB + 2 * beam.dof, 1);
x0(1) = 1;
odeOptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, ...
    'MaxStep', outputStep);
tic;
[tVec, xHist] = ode15s( ...
    @(t, x) hinged_floating_rhs( ...
        t, x, model, NB, beam.dof, hingeTorque), ...
    tOutput, x0, odeOptions);
elapsed = toc;

quaternionHist = xHist(:, 1:4);
quaternionHist = quaternionHist ./ vecnorm(quaternionHist, 2, 2);
basePosition = xHist(:, 5:7);
qHist = xHist(:, 7 + (1:NB));
etaHist = xHist(:, 7 + NB + (1:beam.dof));
hingeAngle = qHist(:, 2);
beamTipElastic = etaHist * beamTipDisplacement.';

nTime = numel(tVec);
baseTheta = zeros(nTime, 1);
beamTipGlobal = zeros(nTime, 3);
for k = 1:nTime
    A = get_rotation_matrix_quaternion(quaternionHist(k, :).');
    Ah = get_rotation_matrix([0; 0; 1], hingeAngle(k));
    baseTheta(k) = atan2(A(2, 1), A(1, 1));
    localTip = [beam.length; 0; 0] + beamTipElastic(k, :).';
    beamTipGlobal(k, :) = (basePosition(k, :).' + A * Ah * localTip).';
end

outputCsv = fullfile(dataDir, 'matlab_hinged_floating_beam.csv');
writematrix_with_header(outputCsv, ...
    [tVec, hingeAngle, beamTipElastic(:, 1:2), basePosition(:, 1:2), ...
     baseTheta, beamTipGlobal(:, 1:2)], ...
    ['time,hinge_angle,beam_tip_elastic_x,beam_tip_elastic_y,', ...
     'base_x,base_y,base_theta,beam_tip_global_x,beam_tip_global_y']);

%% Three-solver comparison
exudynCsv = fullfile(dataDir, 'exudyn_hinged_floating_beam.csv');
chronoCsv = fullfile(dataDir, 'chrono_hinged_floating_beam.csv');
if isfile(exudynCsv) && isfile(chronoCsv)
    exudynData = readmatrix(exudynCsv, 'NumHeaderLines', 1);
    chronoData = readmatrix(chronoCsv, 'NumHeaderLines', 1);
    exudynInterp = interp1(exudynData(:, 1), exudynData(:, 2:9), ...
        tVec, 'linear', 'extrap');
    chronoInterp = interp1(chronoData(:, 1), chronoData(:, 2:9), ...
        tVec, 'linear', 'extrap');

    fig = figure('Color', 'w', 'Position', [100, 100, 1100, 760]);
    layout = tiledlayout(fig, 2, 2, 'Padding', 'compact', ...
        'TileSpacing', 'compact');
    colors = [0.00, 0.45, 0.74; 0.85, 0.33, 0.10; 0.47, 0.67, 0.19];
    plot_comparison(nexttile(layout), tVec, hingeAngle, ...
        exudynInterp(:, 1), chronoInterp(:, 1), ...
        'Hinge angle (rad)', colors, true);
    plot_comparison(nexttile(layout), tVec, basePosition(:, 1), ...
        exudynInterp(:, 4), chronoInterp(:, 4), ...
        'Base x (m)', colors, false);
    plot_comparison(nexttile(layout), tVec, baseTheta, ...
        exudynInterp(:, 6), chronoInterp(:, 6), ...
        'Base angle (rad)', colors, false);
    plot_comparison(nexttile(layout), tVec, beamTipGlobal(:, 2), ...
        exudynInterp(:, 8), chronoInterp(:, 8), ...
        'Beam-tip global y (m)', colors, false);
    exportgraphics(fig, ...
        fullfile(dataDir, 'hinged_floating_beam_comparison.png'), ...
        'Resolution', 220);
else
    warning(['Reference data are missing. Run ', ...
        'refs/exudyn_hinged_floating_beam.py and ', ...
        'refs/chrono_hinged_floating_beam.py, then rerun this script.']);
end

fprintf('MATLAB hinged floating-beam simulation completed in %.3f s.\n', elapsed);
fprintf('Data saved to: %s\n', outputCsv);

function dx = hinged_floating_rhs(~, x, model, NB, dof, hingeTorque)
positionSize = 7 + NB + dof;
theta = x(1:4);
theta = theta / norm(theta);
q = x(7 + (1:NB));
eta = x(7 + NB + (1:dof));
baseVelocity = x(positionSize + (1:6));
qDot = x(positionSize + 6 + (1:NB));
etaDot = x(positionSize + 6 + NB + (1:dof));

A = get_rotation_matrix_quaternion(theta);
baseState.velocity = baseVelocity;
baseState.rotation = A;
tau = [0; hingeTorque];
[qDDot, etaDDot, baseAcceleration] = FD_floating( ...
    model, q, qDot, tau, eta, etaDot, [], baseState);

[~, Gbar] = get_quaternion_G_matrices(theta);
thetaDot = 0.25 * Gbar.' * baseVelocity(1:3);
positionDot = A * baseVelocity(4:6);
dx = [thetaDot; positionDot; qDot; etaDot; ...
      baseAcceleration; qDDot; etaDDot{2}(:)];
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
