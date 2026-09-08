%% Floating base, rigidly attached flexible beam, and a hinged rigid rod
clear; close all; clc;

thisFile = mfilename('fullpath');
testDir = fileparts(thisFile);
srcDir = fileparts(fileparts(testDir));
projectDir = fileparts(srcDir);
dataDir = fullfile(projectDir, 'data', 'base_beam_hinged_rod');
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

addpath(fullfile(srcDir, 'dynamics'), fullfile(srcDir, 'classes'), ...
    fullfile(srcDir, 'utils'), fullfile(testDir, 'support'));

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
beam.geometric_stiffening = true;
beam.geometric_stiffness_plane = 'xy';

femBeam = FEM_Beam_GeometricStiffness(beam);
[~, K3d] = femBeam.get_matrices_3d();
freeDofs = femBeam.get_free_dofs();

% Retain planar [ux, uy, rz] coordinates at every free node.
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
[u0, Sbeam, weights] = femBeam.build_ffrf_data(shapeMatrix, 4, 3, 3);

beamParams.density = beam.density;
beamParams.u0 = u0;
beamParams.S = Sbeam;
beamParams.weights = weights;
beamParams.Kff = Kff;
beamParams.Dff = Dff;
beamParams.geometric_stiffening = beam.geometric_stiffening;
if beam.geometric_stiffening
    geometricData = femBeam.get_geometric_stiffness_data(shapeMatrix);
    beamParams.geometric_C = geometricData.geometric_C;
    beamParams.geometric_D = geometricData.geometric_D;
    beamParams.geometric_rotation_axis = ...
        geometricData.geometric_rotation_axis;
    % Equation (56) uses the prescribed base acceleration along the local
    % beam axis.  It is zero here; rotational stiffening omega^2*D remains.
    beamParams.axial_base_acceleration = 0;
end
beamParams.q = [1; 0; 0; 0; zeros(3 + beam.dof, 1)];
beamParams.qdot = zeros(7 + beam.dof, 1);
flexibleBeam = FlexBody(beamParams);

tipNodeDofs = femBeam.get_node_dofs(beam.elements_num + 1);
tipShape6d = shapeMatrix(tipNodeDofs, :);
beamTipInterface = [tipShape6d(4:6, :); tipShape6d(1:3, :)];
beamTipDisplacement = femBeam.evaluate_displacement( ...
    shapeMatrix, beam.length, 0, 0);

%% Floating rigid base
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

%% Rigid rod hinged at the beam tip
rod.length = 0.20;
rod.mass = 1.5093;
rod.com = [0.10; 0; 0];
rod.Izz = 0.015155;
rodParams.mass = rod.mass;
rodParams.com = rod.com;
% Nonzero transverse inertias keep the spatial floating solve regular.
rodParams.inertia = [rod.Izz; rod.Izz; rod.Izz];
rodBody = RigidBody(rodParams);

%% Flexible-rigid tree
NB = 3;
model.NB = NB;
model.parent = [0; 1; 2];
model.jtype = {'Rz'; 'Rz'; 'Rz'};
% Link 1 is the floating-base placeholder and link 2 is fixed to the base.
model.fixed_joints = [1, 2];
model.Xtree = {eye(6); eye(6); xlt([beam.length, 0, 0])};
model.flexbody = {[]; flexibleBeam; []};
model.rigidbody = {baseBody; []; rodBody};
model.link_type = {'rigid'; 'flexible'; 'rigid'};
model.Phi = {zeros(6, 0); beamTipInterface; zeros(6, 0)};
model.gravity = zeros(6, 1);
model.floating_base = true;

%% Constant internal hinge torque and integration
hingeTorque = 0.02;     % N*m
simulationTime = 5.0;   % s
outputStep = 1.0e-3;    % s
tOutput = (0:outputStep:simulationTime).';

% x = [base quaternion(4); base position(3); q; eta;
%      base spatial velocity(6); qDot; etaDot].
positionSize = 7 + NB + beam.dof;
x0 = zeros(positionSize + 6 + NB + beam.dof, 1);
x0(1) = 1;
odeOptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, ...
    'MaxStep', outputStep);

tic;
[tVec, xHist] = ode15s( ...
    @(t, x) base_beam_hinged_rod_rhs( ...
        t, x, model, NB, beam.dof, hingeTorque), ...
    tOutput, x0, odeOptions);
elapsed = toc;

quaternionHist = xHist(:, 1:4);
quaternionHist = quaternionHist ./ vecnorm(quaternionHist, 2, 2);
basePosition = xHist(:, 5:7);
qHist = xHist(:, 7 + (1:NB));
etaHist = xHist(:, 7 + NB + (1:beam.dof));
hingeAngle = unwrap(qHist(:, 3));
beamTipElastic = etaHist * beamTipDisplacement.';

nTime = numel(tVec);
baseTheta = zeros(nTime, 1);
beamTipGlobal = zeros(nTime, 3);
rodTipGlobal = zeros(nTime, 3);
for k = 1:nTime
    A = get_rotation_matrix_quaternion(quaternionHist(k, :).');
    tipDelta = beamTipInterface * etaHist(k, :).';
    beamEndRotation = get_rotation_matrix([0; 0; 1], tipDelta(3));
    jointRotation = get_rotation_matrix([0; 0; 1], hingeAngle(k));

    baseTheta(k) = atan2(A(2, 1), A(1, 1));
    localBeamTip = [beam.length; 0; 0] + tipDelta(4:6);
    beamTipGlobal(k, :) = (basePosition(k, :).' + A * localBeamTip).';
    rodTipGlobal(k, :) = (beamTipGlobal(k, :).' + ...
        A * beamEndRotation * jointRotation * [rod.length; 0; 0]).';
end
baseTheta = unwrap(baseTheta);

outputCsv = fullfile(dataDir, 'matlab_base_beam_hinged_rod.csv');
writematrix_with_header(outputCsv, ...
    [tVec, hingeAngle, baseTheta, beamTipGlobal(:, 2), ...
     rodTipGlobal(:, 2), beamTipElastic(:, 2), basePosition(:, 1:2), ...
     beamTipGlobal(:, 1), rodTipGlobal(:, 1)], ...
    ['time,hinge_angle,base_theta,beam_tip_global_y,rod_tip_global_y,', ...
     'beam_tip_elastic_y,base_x,base_y,beam_tip_global_x,rod_tip_global_x']);

%% MATLAB, Exudyn, and Chrono comparison
exudynCsv = fullfile(dataDir, 'exudyn_base_beam_hinged_rod.csv');
chronoCsv = fullfile(dataDir, 'chrono_base_beam_hinged_rod.csv');
if isfile(exudynCsv) && isfile(chronoCsv)
    exudynData = read_reference_data(exudynCsv, tVec);
    chronoData = read_reference_data(chronoCsv, tVec);
    matlabData = [hingeAngle, baseTheta, beamTipGlobal(:, 2), ...
        rodTipGlobal(:, 2)];

    comparisonFig = figure('Color', 'w', 'Position', [100, 100, 1100, 760]);
    layout = tiledlayout(comparisonFig, 2, 2, ...
        'Padding', 'compact', 'TileSpacing', 'compact');
    colors = [0.00, 0.45, 0.74; ...
              0.85, 0.33, 0.10; ...
              0.47, 0.67, 0.19];
    labels = {'Hinge angle (rad)', 'Base angle (rad)', ...
        'Beam-tip global y (m)', 'Rod-tip global y (m)'};
    for i = 1:4
        plot_comparison(nexttile(layout), tVec, matlabData(:, i), ...
            exudynData(:, i), chronoData(:, i), labels{i}, colors, i == 1);
    end

    comparisonPng = fullfile(dataDir, ...
        'base_beam_hinged_rod_comparison.png');
    exportgraphics(comparisonFig, comparisonPng, 'Resolution', 220);

    exudynRms = sqrt(mean((matlabData - exudynData).^2, 1));
    chronoRms = sqrt(mean((matlabData - chronoData).^2, 1));
    fprintf(['Exudyn RMS [hinge, base, beam tip y, rod tip y]: ', ...
        '[%.4e rad, %.4e rad, %.4e m, %.4e m]\n'], exudynRms);
    fprintf(['Chrono RMS [hinge, base, beam tip y, rod tip y]: ', ...
        '[%.4e rad, %.4e rad, %.4e m, %.4e m]\n'], chronoRms);
    fprintf('Comparison figure: %s\n', comparisonPng);
else
    warning(['Reference data are missing. Run ', ...
        'refs/exudyn_base_beam_hinged_rod.py and ', ...
        'refs/chrono_base_beam_hinged_rod.py, then rerun this script.']);
end

fprintf('MATLAB base-beam-hinged-rod simulation completed in %.3f s.\n', ...
    elapsed);
fprintf('Data saved to: %s\n', outputCsv);

function dx = base_beam_hinged_rod_rhs( ...
        ~, x, model, NB, dof, hingeTorque)
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
tau = [0; 0; hingeTorque];
[qDDot, etaDDot, baseAcceleration] = FD_floating( ...
    model, q, qDot, tau, eta, etaDot, [], baseState);

[~, Gbar] = get_quaternion_G_matrices(theta);
thetaDot = 0.25 * Gbar.' * baseVelocity(1:3);
positionDot = A * baseVelocity(4:6);
dx = [thetaDot; positionDot; qDot; etaDot; ...
      baseAcceleration; qDDot; etaDDot{2}(:)];
end

function referenceInterp = read_reference_data(csvFile, tVec)
data = readmatrix(csvFile, 'NumHeaderLines', 1);
data = data(all(isfinite(data), 2), :);
data = sortrows(data, 1);
[time, uniqueIds] = unique(data(:, 1), 'stable');
data = data(uniqueIds, :);
referenceInterp = interp1(time, data(:, 2:5), tVec, 'linear', 'extrap');
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
