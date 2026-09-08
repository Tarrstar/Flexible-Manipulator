%% Fixed-base flexible beam manipulator
clear; close all; clc;

thisFile = mfilename('fullpath');
testDir = fileparts(thisFile);
srcDir = fileparts(fileparts(testDir));
projectDir = fileparts(srcDir);
dataDir = fullfile(projectDir, 'data', 'fixed_FBM');
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
[M3d, K3d] = femBeam.get_matrices_3d();
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

%% Four rigid links and three actuated revolute joints
rigidMass = [1.5093, 0.1029, 0.1139, 0.1009];
rigidComX = [0.1, 0.0709, 0.0723, 0.0504];
rigidIzz = [0.015155, 0.002069, 0.002387, 0.001107];
rigidBodies = cell(4, 1);
for i = 1:4
    p.mass = rigidMass(i);
    p.com = [rigidComX(i); 0; 0];
    p.inertia = [rigidIzz(i); rigidIzz(i); rigidIzz(i)];
    rigidBodies{i} = RigidBody(p);
end

NB = 5;
model.NB = NB;
model.parent = [0; 1; 2; 3; 4];
model.jtype = {'Rz'; 'Rz'; 'Rz'; 'Rz'; 'Rz'};
% The beam root and beam-tip/first-link connection are rigidly fixed.
model.fixed_joints = [1, 2];
model.Xtree = { ...
    eye(6); ...
    xlt([beam.length, 0, 0]); ...
    xlt([0.2, 0, 0]); ...
    xlt([0.1418, 0, 0]); ...
    xlt([0.1448, 0, 0])};
model.flexbody = [{flexibleBeam}; repmat({[]}, 4, 1)];
model.rigidbody = [{[]}; rigidBodies];
model.link_type = [{'flexible'}; repmat({'rigid'}, 4, 1)];
model.Phi = [{beamTipInterface}; repmat({zeros(6, 0)}, 4, 1)];
model.gravity = zeros(6, 1);

%% Integration
simulationTime = 5.0;
outputStep = 1.0e-3;
jointTorque = 0.02;
tOutput = (0:outputStep:simulationTime).';
tau = [0; 0; jointTorque; jointTorque; jointTorque];

% x = [q; eta; qDot; etaDot].
x0 = zeros(2 * (NB + beam.dof), 1);
odeOptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, ...
    'MaxStep', outputStep);
tic;
[tVec, xHist] = ode15s( ...
    @(t, x) fixed_base_rhs(t, x, model, NB, beam.dof, 1, tau), ...
    tOutput, x0, odeOptions);
elapsed = toc;

qHist = xHist(:, 1:NB);
jointAngles = qHist(:, 3:5);
etaHist = xHist(:, NB + (1:beam.dof));
beamTipElastic = etaHist * beamTipDisplacement.';
beamTipGlobal = beamTipElastic + repmat([beam.length, 0, 0], numel(tVec), 1);

outputCsv = fullfile(dataDir, 'matlab_fixed_FBM.csv');
writematrix_with_header(outputCsv, ...
    [tVec, jointAngles, beamTipElastic(:, 1:2), beamTipGlobal(:, 1:2)], ...
    ['time,joint_1_angle,joint_2_angle,joint_3_angle,', ...
     'beam_tip_elastic_x,beam_tip_elastic_y,', ...
     'beam_tip_global_x,beam_tip_global_y']);

%% Three-solver comparison
exudynCsv = fullfile(dataDir, 'exudyn_fixed_FBM.csv');
chronoCsv = fullfile(dataDir, 'chrono_fixed_FBM.csv');
if isfile(exudynCsv) && isfile(chronoCsv)
    exudynData = readmatrix(exudynCsv, 'NumHeaderLines', 1);
    chronoData = readmatrix(chronoCsv, 'NumHeaderLines', 1);
    exudynJoint = interp1(exudynData(:, 1), exudynData(:, 2:4), ...
        tVec, 'linear', 'extrap');
    chronoJoint = interp1(chronoData(:, 1), chronoData(:, 2:4), ...
        tVec, 'linear', 'extrap');
    exudynTipY = interp1(exudynData(:, 1), exudynData(:, 6), ...
        tVec, 'linear', 'extrap');
    chronoTipY = interp1(chronoData(:, 1), chronoData(:, 6), ...
        tVec, 'linear', 'extrap');

    fig = figure('Color', 'w', 'Position', [100, 100, 1100, 760]);
    layout = tiledlayout(fig, 2, 2, 'Padding', 'compact', ...
        'TileSpacing', 'compact');
    colors = [0.00, 0.45, 0.74; 0.85, 0.33, 0.10; 0.47, 0.67, 0.19];
    for i = 1:3
        ax = nexttile(layout);
        plot(ax, tVec, jointAngles(:, i), '-', 'Color', colors(1, :), ...
            'LineWidth', 1.4); hold(ax, 'on');
        plot(ax, tVec, exudynJoint(:, i), '--', 'Color', colors(2, :), ...
            'LineWidth', 1.3);
        plot(ax, tVec, chronoJoint(:, i), ':', 'Color', colors(3, :), ...
            'LineWidth', 1.6); hold(ax, 'off');
        grid(ax, 'on'); box(ax, 'on');
        xlabel(ax, 'Time (s)'); ylabel(ax, sprintf('Joint %d (rad)', i));
        if i == 1
            legend(ax, {'MATLAB FFRF', 'Exudyn', 'Chrono'}, ...
                'Location', 'best');
        end
    end
    ax = nexttile(layout);
    plot(ax, tVec, beamTipElastic(:, 2), '-', 'Color', colors(1, :), ...
        'LineWidth', 1.4); hold(ax, 'on');
    plot(ax, tVec, exudynTipY, '--', 'Color', colors(2, :), ...
        'LineWidth', 1.3);
    plot(ax, tVec, chronoTipY, ':', 'Color', colors(3, :), ...
        'LineWidth', 1.6); hold(ax, 'off');
    grid(ax, 'on'); box(ax, 'on');
    xlabel(ax, 'Time (s)'); ylabel(ax, 'Beam-tip elastic y (m)');
    exportgraphics(fig, fullfile(dataDir, 'fixed_FBM_comparison.png'), ...
        'Resolution', 220);
else
    warning(['Reference data are missing. Run refs/exudyn_fixed_FBM.py ', ...
        'and refs/chrono_fixed_FBM.py, then rerun this script.']);
end

fprintf('MATLAB fixed-FBM simulation completed in %.3f s.\n', elapsed);
fprintf('Data saved to: %s\n', outputCsv);
