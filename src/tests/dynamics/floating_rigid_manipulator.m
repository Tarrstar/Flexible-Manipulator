%% Floating rigid-base three-joint robot: FD_floating vs spatial_v2 FDfb
clear; close all; clc;

thisFile = mfilename('fullpath');
testDir = fileparts(thisFile);
srcDir = fileparts(fileparts(testDir));
projectDir = fileparts(srcDir);
dataDir = fullfile(projectDir, 'data', 'floating_rigid_manipulator');
if ~exist(dataDir, 'dir'), mkdir(dataDir); end

spatialDir = fullfile(fileparts(projectDir), 'spatial_v2.1');
if ~isfile(fullfile(spatialDir, 'dynamics', 'FDfb.m'))
    error('spatial_v2.1 was not found at: %s', spatialDir);
end

addpath(genpath(spatialDir), fullfile(srcDir, 'classes'), ...
    fullfile(srcDir, 'dynamics'), fullfile(srcDir, 'utils'));

%% Equivalent rigid floating base
beam.length = 1.0;
beam.width = 0.032;
beam.thickness = 0.0015;
beam.density = 7.542e3;

beamMass = beam.density * beam.length * beam.width * beam.thickness;
beamCom = [0.5 * beam.length; 0; 0];
beamInertia = diag(beamMass / 12 * [ ...
    beam.width^2 + beam.thickness^2; ...
    beam.length^2 + beam.thickness^2; ...
    beam.length^2 + beam.width^2]);

fixedLength = 0.2;
fixedMass = 1.5093;
fixedCom = [beam.length + 0.1; 0; 0];
fixedInertia = diag([0.015155, 0.015155, 0.015155]);

baseMass = beamMass + fixedMass;
baseCom = (beamMass * beamCom + fixedMass * fixedCom) / baseMass;
beamComSkew = get_skew_sym_matrix(beamCom);
fixedComSkew = get_skew_sym_matrix(fixedCom);
baseComSkew = get_skew_sym_matrix(baseCom);
baseInertiaOrigin = beamInertia - beamMass * beamComSkew^2 + ...
    fixedInertia - fixedMass * fixedComSkew^2;
baseInertia = baseInertiaOrigin + baseMass * baseComSkew^2;
baseInertia = 0.5 * (baseInertia + baseInertia.');

params.mass = baseMass;
params.com = baseCom;
params.inertia = baseInertia;
rigidBodies = cell(4, 1);
rigidBodies{1} = RigidBody(params);

linkMass = [0.1029, 0.1139, 0.1009];
linkComX = [0.0709, 0.0723, 0.0504];
linkIzz = [0.002069, 0.002387, 0.001107];
for i = 1:3
    params.mass = linkMass(i);
    params.com = [linkComX(i); 0; 0];
    params.inertia = linkIzz(i) * eye(3);
    rigidBodies{i + 1} = RigidBody(params);
end

%% FD_floating model: body 1 is the floating base
NB = 4;
model.NB = NB;
model.parent = [0; 1; 2; 3];
model.jtype = {'Rz'; 'Rz'; 'Rz'; 'Rz'};
model.fixed_joints = 1;
model.Xtree = { ...
    eye(6); ...
    xlt([beam.length + fixedLength, 0, 0]); ...
    xlt([0.1418, 0, 0]); ...
    xlt([0.1448, 0, 0])};
model.rigidbody = rigidBodies;
model.flexbody = repmat({[]}, NB, 1);
model.link_type = repmat({'rigid'}, NB, 1);
model.Phi = repmat({zeros(6, 0)}, NB, 1);
model.gravity = zeros(6, 1);
model.floating_base = true;

%% spatial_v2 model constructed from exactly the same inertias
spatialSeed.NB = NB;
spatialSeed.parent = model.parent.';
spatialSeed.jtype = model.jtype.';
spatialSeed.Xtree = model.Xtree.';
spatialSeed.I = cell(1, NB);
for i = 1:NB
    spatialSeed.I{i} = rigidBodies{i}.get_spatial_inertia();
end
spatialSeed.gravity = [0, 0, 0];
spatialModel = floatbase(spatialSeed);

%% Time integration
simulationTime = 0.5;
outputStep = 0.002;
tOutput = (0:outputStep:simulationTime).';
jointTorque = 0.1 * ones(3, 1);
odeOptions = odeset('RelTol', 1e-9, 'AbsTol', 1e-11, ...
    'InitialStep', 1e-7, 'MaxStep', outputStep);

% FD_floating state: [quaternion; position; body velocity; q; qdot].
x0 = zeros(13 + 2 * NB, 1);
x0(1) = 1;
[tFd, xFd] = ode45( ...
    @(t, x) fd_floating_rhs(t, x, model, NB, jointTorque), ...
    tOutput, x0, odeOptions);

% spatial_v2 state: [xfb(13); q(3); qdot(3)].
x0Spatial = zeros(13 + 2 * numel(jointTorque), 1);
x0Spatial(1) = 1;
[tSpatial, xSpatial] = ode45( ...
    @(t, x) spatial_v2_rhs(t, x, spatialModel, jointTorque), ...
    tOutput, x0Spatial, odeOptions);

%% Common output coordinates
fdQuaternion = normalize_quaternion_rows(xFd(:, 1:4));
spatialQuaternion = normalize_quaternion_rows(xSpatial(:, 1:4));
fdPosition = xFd(:, 5:7);
spatialPosition = xSpatial(:, 5:7);
fdJoint = xFd(:, 13 + (2:NB));
spatialJoint = xSpatial(:, 13 + (1:3));

fdBaseAngle = zeros(numel(tFd), 1);
spatialBaseAngle = zeros(numel(tSpatial), 1);
for k = 1:numel(tFd)
    A = get_rotation_matrix_quaternion(fdQuaternion(k, :).');
    fdBaseAngle(k) = atan2(A(2, 1), A(1, 1));

    % spatial_v2 quaternion maps inertial coordinates into the body frame.
    E = rq(spatialQuaternion(k, :).');
    ASpatial = E.';
    spatialBaseAngle(k) = atan2(ASpatial(2, 1), ASpatial(1, 1));
end
fdBaseAngle = unwrap(fdBaseAngle);
spatialBaseAngle = unwrap(spatialBaseAngle);

basePositionError = fdPosition - spatialPosition;
baseAngleError = fdBaseAngle - spatialBaseAngle;
jointError = fdJoint - spatialJoint;

outputCsv = fullfile(dataDir, 'fd_floating_vs_spatial_v2.csv');
writematrix_with_header(outputCsv, ...
    [tFd, fdPosition(:, 1:2), spatialPosition(:, 1:2), ...
     fdBaseAngle, spatialBaseAngle, fdJoint, spatialJoint], ...
    ['time,fd_base_x,fd_base_y,spatial_base_x,spatial_base_y,', ...
     'fd_base_theta,spatial_base_theta,', ...
     'fd_joint_1,fd_joint_2,fd_joint_3,', ...
     'spatial_joint_1,spatial_joint_2,spatial_joint_3']);

%% Comparison figure
comparisonFig = figure('Color', 'w', 'Position', [80, 80, 1150, 850]);
layout = tiledlayout(comparisonFig, 3, 2, ...
    'Padding', 'compact', 'TileSpacing', 'compact');
colors = [0.00, 0.45, 0.74; 0.85, 0.33, 0.10];

plot_pair(nexttile(layout), tFd, fdPosition(:, 1), ...
    tSpatial, spatialPosition(:, 1), 'Base x (m)', colors, true);
plot_pair(nexttile(layout), tFd, fdPosition(:, 2), ...
    tSpatial, spatialPosition(:, 2), 'Base y (m)', colors, false);
plot_pair(nexttile(layout), tFd, fdBaseAngle, ...
    tSpatial, spatialBaseAngle, 'Base angle (rad)', colors, false);
for i = 1:3
    plot_pair(nexttile(layout), tFd, fdJoint(:, i), ...
        tSpatial, spatialJoint(:, i), sprintf('Joint %d (rad)', i), ...
        colors, false);
end

comparisonPng = fullfile(dataDir, 'fd_floating_vs_spatial_v2.png');
exportgraphics(comparisonFig, comparisonPng, 'Resolution', 220);

fprintf('Floating rigid robot: FD_floating vs spatial_v2 FDfb\n');
fprintf('  max base-position error = %.6e m\n', ...
    max(abs(basePositionError), [], 'all'));
fprintf('  max base-angle error    = %.6e rad\n', max(abs(baseAngleError)));
fprintf('  max joint-angle errors  = [%s] rad\n', ...
    sprintf('%.6e ', max(abs(jointError), [], 1)));
fprintf('  data: %s\n', outputCsv);
fprintf('  figure: %s\n', comparisonPng);

function dx = fd_floating_rhs(~, x, model, NB, jointTorque)
theta = x(1:4);
theta = theta / norm(theta);
baseVelocity = x(8:13);
q = x(13 + (1:NB));
qd = x(13 + NB + (1:NB));

A = get_rotation_matrix_quaternion(theta);
baseState.velocity = baseVelocity;
baseState.rotation = A;
tau = [0; jointTorque];
[qdd, ~, baseAccel] = FD_floating( ...
    model, q, qd, tau, [], [], [], baseState);

[~, Gbar] = get_quaternion_G_matrices(theta);
thetaDot = 0.25 * Gbar.' * baseVelocity(1:3);
positionDot = A * baseVelocity(4:6);
baseVelocityDot = baseAccel;

dx = [thetaDot; positionDot; baseVelocityDot; qd; qdd];
end

function dx = spatial_v2_rhs(~, x, model, jointTorque)
nJoint = numel(jointTorque);
xfb = x(1:13);
xfb(1:4) = xfb(1:4) / norm(xfb(1:4));
q = x(13 + (1:nJoint));
qd = x(13 + nJoint + (1:nJoint));
[xfbDot, qdd] = FDfb(model, xfb, q, qd, jointTorque);
dx = [xfbDot; qd; qdd];
end

function quaternion = normalize_quaternion_rows(quaternion)
quaternion = quaternion ./ vecnorm(quaternion, 2, 2);
end

function plot_pair(ax, t1, y1, t2, y2, yLabelText, colors, showLegend)
plot(ax, t1, y1, '-', 'Color', colors(1, :), 'LineWidth', 1.4);
hold(ax, 'on');
plot(ax, t2, y2, '--', 'Color', colors(2, :), 'LineWidth', 1.3);
hold(ax, 'off');
grid(ax, 'on'); box(ax, 'on');
xlabel(ax, 'Time (s)'); ylabel(ax, yLabelText);
if showLegend
    legend(ax, {'FD\_floating', 'spatial\_v2 FDfb'}, 'Location', 'best');
end
end
