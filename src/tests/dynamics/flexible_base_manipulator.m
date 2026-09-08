%% Floating rigid base, flexible beam, and three-joint manipulator
clear; close all; clc;

thisFile = mfilename('fullpath');
testDir = fileparts(thisFile);
srcDir = fileparts(fileparts(testDir));
projectDir = fileparts(srcDir);
dataDir = fullfile(projectDir, 'data', 'flexible_base_manipulator');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

addpath(fullfile(srcDir, 'dynamics'), fullfile(srcDir, 'classes'), ...
    fullfile(srcDir, 'utils'), fullfile(srcDir, 'integrators'), ...
    fullfile(testDir, 'support'));

%% Beam
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
[M3d, K3d] = femBeam.get_matrices_3d();
freeDofs = femBeam.get_free_dofs();

% Retain linear axial displacement and planar bending coordinates
% [ux, uy, rz] at every free node. Kff remains the linear FE stiffness;
% no geometric stiffness or nonlinear axial-bending strain is included.
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
Mff = full(M3d(beamCoordIds, beamCoordIds));
Dff = 2 * beam.zeta / beam.omega_n * Kff;

% Each element has 4*3*3 nodes for Gauss quadrature, the total number of
% the nodes is 4*3*3*5.
[u0, Sbeam, weights] = femBeam.build_ffrf_data(shapeMatrix, 4, 3, 3);
beamParams.density = beam.density;
beamParams.u0 = u0;
beamParams.S = Sbeam;
beamParams.weights = weights;
% beamParams.Mff = Mff;
beamParams.Kff = Kff;
beamParams.Dff = Dff;
beamParams.q = [1; 0; 0; 0; zeros(3 + beam.dof, 1)];
beamParams.qdot = zeros(7 + beam.dof, 1);
flexibleBase = FlexBody(beamParams);

tipNodeDofs = femBeam.get_node_dofs(beam.elements_num + 1);
tipShape6d = shapeMatrix(tipNodeDofs, :);
beamTipInterface = [tipShape6d(4:6, :); tipShape6d(1:3, :)];
beamTipDisplacement = femBeam.evaluate_displacement( ...
    shapeMatrix, beam.length, 0, 0);

%% Floating base
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

%% Rigid manipulator
rigidMass = [1.5093, 0.1029, 0.1139, 0.1009];
rigidComX = [0.1, 0.0709, 0.0723, 0.0504];
rigidIzz = [0.015155, 0.002069, 0.002387, 0.001107];
rigidBodies = cell(4, 1);
for i = 1:4
    p.mass = rigidMass(i);
    p.com = [rigidComX(i); 0; 0];
    % FD_flex uses a 3-D spatial inertia even for this planar example.
    % Zero transverse inertias make a floating base numerically singular.
    p.inertia = [rigidIzz(i); rigidIzz(i); rigidIzz(i)];
    rigidBodies{i} = RigidBody(p);
end

%% Flexible-Rigid coupling system
NB = 6;
model.NB = NB;
model.parent = [0; 1; 2; 3; 4; 5];
model.jtype = {'Rz'; 'Rz'; 'Rz'; 'Rz'; 'Rz'; 'Rz'};
model.fixed_joints = [1, 2, 3];
model.Xtree = { ...
    eye(6); ...
    eye(6); ...
    xlt([beam.length, 0, 0]); ...
    xlt([0.2, 0, 0]); ...
    xlt([0.1418, 0, 0]); ...
    xlt([0.1448, 0, 0])};
model.flexbody = [{[]}; {flexibleBase}; repmat({[]}, 4, 1)];
model.rigidbody = [{baseBody}; {[]}; rigidBodies];
model.link_type = [{'rigid'}; {'flexible'}; repmat({'rigid'}, 4, 1)];
model.Phi = [{zeros(6, 0)}; {beamTipInterface}; ...
    repmat({zeros(6, 0)}, 4, 1)];
model.gravity = zeros(6, 1);
model.floating_base = true;

%% Implicit integration
simulationTime = 5.0;
outputStep = 0.001;
jointTorqueAmplitude = 0.02;
jointTorqueFrequency = 1.0;
jointTorque = @(t) jointTorqueAmplitude * ...
    sin(2 * pi * jointTorqueFrequency * t);
nOutputSteps = round(simulationTime / outputStep);
tOutput = linspace(0, simulationTime, nOutputSteps + 1).';

% x = [base quaternion(4); base position(3); q; eta;
%      base spatial velocity(6); qDot; etaDot].
nState = 13 + 2 * NB + 2 * beam.dof;
x0 = zeros(nState, 1);
x0(1) = 1;
% Original generalized-alpha integration:
% integratorOptions.RelTol = 1e-6;
% integratorOptions.AbsTol = 1e-8;
% integratorOptions.MaxStep = 1e-3;
% integratorOptions.SpectralRadius = 0.9;
% integratorOptions.PositionSize = 7 + NB + beam.dof;
% integratorOptions.QuaternionIndices = 1:4;
% [tVec, xHist] = generalized_alpha_integrator( ...
%     @(t, x) floating_base_rhs( ...
%         t, x, model, NB, beam.dof, jointTorque), ...
%     tOutput, x0, integratorOptions);

odeOptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, ...
    'MaxStep', outputStep);
tic;
lastwarn('');
[tVec, xHist] = ode15s( ...
    @(t, x) floating_base_rhs( ...
        t, x, model, NB, beam.dof, jointTorque), ...
    tOutput, x0, odeOptions);
toc;

quaternionHist = xHist(:, 1:4);
quaternionHist = quaternionHist ./ vecnorm(quaternionHist, 2, 2);
basePosition = xHist(:, 5:7);
qHist = xHist(:, 7 + (1:NB));
etaHist = xHist(:, 7 + NB + (1:beam.dof));
jointAngles = qHist(:, 4:6);
beamTipElastic = etaHist * beamTipDisplacement.';
beamTipElasticX = beamTipElastic(:, 1);
beamTipElasticY = beamTipElastic(:, 2);

nTime = numel(tVec);
baseTheta = zeros(nTime, 1);
beamTipGlobal = zeros(nTime, 3);
for k = 1:nTime
    A = get_rotation_matrix_quaternion(quaternionHist(k, :).');
    baseTheta(k) = atan2(A(2, 1), A(1, 1));
    localTip = [beam.length; 0; 0] + ...
        beamTipDisplacement * etaHist(k, :).';
    beamTipGlobal(k, :) = (basePosition(k, :).' + A * localTip).';
end

%%
outputCsv = fullfile(dataDir, 'matlab_floating_base.csv');
writematrix_with_header(outputCsv, ...
    [tVec, jointAngles, beamTipElasticY, basePosition(:, 1:2), ...
     baseTheta, beamTipGlobal(:, 1:2), beamTipElasticX], ...
    ['time,joint_1_angle,joint_2_angle,joint_3_angle,', ...
     'beam_tip_elastic_y,base_x,base_y,base_theta,', ...
     'beam_tip_global_x,beam_tip_global_y,beam_tip_elastic_x']);

beamTipFig = figure('Color', 'w', 'Position', [140, 140, 820, 500]);
beamTipAx = axes(beamTipFig);
plot(beamTipAx, tVec, beamTipElasticY, '-', ...
    'Color', [0.00, 0.45, 0.74], 'LineWidth', 1.6);
grid(beamTipAx, 'on');
box(beamTipAx, 'on');
xlabel(beamTipAx, 'Time (s)');
ylabel(beamTipAx, 'Beam-tip elastic y displacement (m)');
title(beamTipAx, 'Beam-tip deformation in the floating frame');

beamTipPng = fullfile(dataDir, 'beam_tip_displacement.png');
exportgraphics(beamTipFig, beamTipPng, 'Resolution', 220);

exudynCsv = fullfile(dataDir, 'exudyn_floating_base.csv');
if ~isfile(exudynCsv)
    error(['Exudyn comparison data is missing. Run refs/', ...
        'exudyn_flexible_base_manipulator.py first.']);
end

exudynData = readmatrix(exudynCsv, 'NumHeaderLines', 1);
exudynData = sortrows(exudynData, 1);
[exudynTime, uniqueIds] = unique(exudynData(:, 1), 'stable');
exudynData = exudynData(uniqueIds, :);
exudynJoint = interp1(exudynTime, exudynData(:, 2:4), ...
    tVec, 'linear', 'extrap');
exudynBaseTheta = interp1(exudynTime, exudynData(:, 8), ...
    tVec, 'linear', 'extrap');

chronoCsv = fullfile(dataDir, 'chrono_floating_base.csv');
if ~isfile(chronoCsv)
    error(['Chrono comparison data is missing. Run refs/', ...
        'chrono_flexible_base_manipulator.py first.']);
end
chronoData = readmatrix(chronoCsv, 'NumHeaderLines', 1);
chronoData = sortrows(chronoData, 1);
[chronoTime, uniqueIds] = unique(chronoData(:, 1), 'stable');
chronoData = chronoData(uniqueIds, :);
chronoJoint = interp1(chronoTime, chronoData(:, 2:4), ...
    tVec, 'linear', 'extrap');
chronoBaseTheta = interp1(chronoTime, chronoData(:, 8), ...
    tVec, 'linear', 'extrap');

comparisonFig = figure('Color', 'w', 'Position', [100, 100, 1100, 760]);
layout = tiledlayout(comparisonFig, 2, 2, ...
    'Padding', 'compact', 'TileSpacing', 'compact');
lineColors = [0.00, 0.45, 0.74; ...
              0.85, 0.33, 0.10; ...
              0.47, 0.67, 0.19];
for i = 1:3
    ax = nexttile(layout);
    plot(ax, tVec, jointAngles(:, i), '-', 'Color', lineColors(1, :), ...
        'LineWidth', 1.4); hold(ax, 'on');
    plot(ax, tVec, exudynJoint(:, i), '--', 'Color', lineColors(2, :), ...
        'LineWidth', 1.3);
    plot(ax, tVec, chronoJoint(:, i), ':', 'Color', lineColors(3, :), ...
        'LineWidth', 1.6); hold(ax, 'off');
    grid(ax, 'on'); box(ax, 'on');
    xlabel(ax, 'Time (s)'); ylabel(ax, sprintf('Joint %d (rad)', i));
    if i == 1
        legend(ax, {'MATLAB FFRF', 'Exudyn', 'Chrono'}, ...
            'Location', 'best');
    end
end

ax = nexttile(layout);
plot(ax, tVec, baseTheta, '-', 'Color', lineColors(1, :), ...
    'LineWidth', 1.4); hold(ax, 'on');
plot(ax, tVec, exudynBaseTheta, '--', 'Color', lineColors(2, :), ...
    'LineWidth', 1.3);
plot(ax, tVec, chronoBaseTheta, ':', 'Color', lineColors(3, :), ...
    'LineWidth', 1.6); hold(ax, 'off');
grid(ax, 'on'); box(ax, 'on');
xlabel(ax, 'Time (s)'); ylabel(ax, 'Floating base angle (rad)');
legend(ax, {'MATLAB FFRF', 'Exudyn', 'Chrono'}, 'Location', 'best');

comparisonPng = fullfile(dataDir, 'floating_base_comparison.png');
exportgraphics(comparisonFig, comparisonPng, 'Resolution', 220);

% jointError = jointAngles - exudynJoint;
% baseThetaError = baseTheta - exudynBaseTheta;
% chronoJointError = jointAngles - chronoJoint;
% chronoBaseThetaError = baseTheta - chronoBaseTheta;
%
% fprintf('Exudyn joint RMS errors: [%s] rad\n', ...
%     sprintf('%.4e ', sqrt(mean(jointError.^2, 1))));
% fprintf('Exudyn floating-base angle RMS error: %.4e rad\n', ...
%     sqrt(mean(baseThetaError.^2)));
% fprintf('Chrono joint RMS errors: [%s] rad\n', ...
%     sprintf('%.4e ', sqrt(mean(chronoJointError.^2, 1))));
% fprintf('Chrono floating-base angle RMS error: %.4e rad\n', ...
%     sqrt(mean(chronoBaseThetaError.^2)));
% fprintf('MATLAB data: %s\n', outputCsv);
% fprintf('Comparison: %s\n', comparisonPng);
