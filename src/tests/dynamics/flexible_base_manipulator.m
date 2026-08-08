%% Flexible-base three-joint manipulator under constant joint torques
clear; close all; clc;

thisFile = mfilename('fullpath');
testDir = fileparts(thisFile);
srcDir = fileparts(fileparts(testDir));
projectDir = fileparts(srcDir);
dataDir = fullfile(projectDir, 'data', 'flexible_base_manipulator');
if ~exist(dataDir, 'dir'), mkdir(dataDir); end

addpath(fullfile(srcDir, 'classes'), fullfile(srcDir, 'dynamics'), ...
    fullfile(srcDir, 'utils'), fullfile(testDir, 'support'));

%% Flexible cantilever: parameters from joint_friction_test.m
beam.length = 1.0;
beam.width = 0.032;
beam.thickness = 0.0015;
beam.density = 7.542e3;
beam.youngs_modulus = 1.983e11;
beam.poisson_ratio = 0.30;
beam.elements_num = 5;
beam.fixed_root = true;
beam.zeta = 0.0053;
beam.omega_n = 7.808;

femBeam = FEM_Beam(beam);
[M3d, K3d] = femBeam.get_matrices_3d();
freeDofs = femBeam.get_free_dofs();

% Match the reference planar Euler-Bernoulli coordinates exactly:
% transverse y displacement and rotation about z at every free node.
nEta = 2 * beam.elements_num;
planarGlobalDofs = zeros(1, nEta);
for node = 2:(beam.elements_num + 1)
    nodeDofs = femBeam.get_node_dofs(node);
    ids = 2 * (node - 2) + (1:2);
    planarGlobalDofs(ids) = nodeDofs([2, 6]);
end
[isFree, planarCoordIds] = ismember(planarGlobalDofs, freeDofs);
if ~all(isFree)
    error('Could not map planar beam DOFs into the fixed-root FE matrices.');
end

shapeMatrix = zeros(femBeam.get_total_dof(), nEta);
for k = 1:nEta
    shapeMatrix(planarGlobalDofs(k), k) = 1;
end

Kff = full(K3d(planarCoordIds, planarCoordIds));
Mff = full(M3d(planarCoordIds, planarCoordIds));
Kff = 0.5 * (Kff + Kff.');
Mff = 0.5 * (Mff + Mff.');
Dff = 2 * beam.zeta * beam.omega_n * Mff;

[u0, Sbeam, weights] = femBeam.build_ffrf_data(shapeMatrix, 4, 3, 3);
rodParams.density = beam.density;
rodParams.u0 = u0;
rodParams.S = Sbeam;
rodParams.weights = weights;
rodParams.Kff = Kff;
rodParams.Dff = Dff;
rodParams.q = [1; 0; 0; 0; zeros(3 + nEta, 1)];
rodParams.qdot = zeros(7 + nEta, 1);
flexibleBase = FlexBody(rodParams);

tipNodeDofs = femBeam.get_node_dofs(beam.elements_num + 1);
tipShape6d = shapeMatrix(tipNodeDofs, :);
beamTipInterface = [tipShape6d(4:6, :); tipShape6d(1:3, :)];
beamTipDisplacement = femBeam.evaluate_displacement( ...
    shapeMatrix, beam.length, 0, 0);

%% Flexible base plus four rigid links and three independent joints
% Reference model generalized coordinates are:
% [beam tip y, beam tip rz, joint 1, joint 2, joint 3].
NB = 5;
model.NB = NB;
model.parent = [0; 1; 2; 3; 4];
model.jtype = {'Rz'; 'Rz'; 'Rz'; 'Rz'; 'Rz'};
model.fixed_joints = [1, 2];
model.Xtree = { ...
    eye(6); ...
    xlt([beam.length, 0, 0]); ...
    xlt([0.2, 0, 0]); ...
    xlt([0.1418, 0, 0]); ...
    xlt([0.1448, 0, 0])};

rigidMass = [1.5093, 0.1029, 0.1139, 0.1009];
rigidComX = [0.1, 0.0709, 0.0723, 0.0504];
rigidIzz = [0.015155, 0.002069, 0.002387, 0.001107];
rigidBodies = cell(4, 1);
for i = 1:4
    p.mass = rigidMass(i);
    p.com = [rigidComX(i); 0; 0];
    p.inertia = [0; 0; rigidIzz(i)];
    rigidBodies{i} = RigidBody(p);
end

model.flexbody = [{flexibleBase}; repmat({[]}, 4, 1)];
model.rigidbody = [{[]}; rigidBodies];
model.link_type = [{'flexible'}; repmat({'rigid'}, 4, 1)];
model.Phi = [{beamTipInterface}; repmat({zeros(6, 0)}, 4, 1)];
model.gravity = zeros(6, 1);

%% ODE45 simulation: 0.1 N m on each independent joint
simulationTime = 5.0;
jointTorque = 0.005;
x0 = zeros(2 * NB + 2 * nEta, 1);
odeOptions = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

[tVec, xHist] = ode45( ...
    @(t, x) flexible_base_rhs(t, x, model, NB, nEta, jointTorque), ...
    [0, simulationTime], x0, odeOptions);

qHist = xHist(:, 1:NB);
etaHist = xHist(:, 2 * NB + (1:nEta));
jointAngles = qHist(:, 3:5);
beamTipY = etaHist * beamTipDisplacement(2, :).';

outputCsv = fullfile(dataDir, 'matlab_joint_and_beam_tip.csv');
writematrix_with_header(outputCsv, ...
    [tVec(:), jointAngles, beamTipY], ...
    'time,joint_1_angle,joint_2_angle,joint_3_angle,beam_tip_y');

%% Joint-angle and beam-tip comparison
figure('Color', 'w', 'Name', 'Flexible-base manipulator response');
subplot(2, 2, 1);
plot(tVec, jointAngles(:, 1), 'LineWidth', 1.3);
grid on;
ylabel('Joint 1 angle (rad)');
title('Joint 1');

subplot(2, 2, 2);
plot(tVec, jointAngles(:, 2), 'LineWidth', 1.3);
grid on;
ylabel('Joint 2 angle (rad)');
title('Joint 2');

subplot(2, 2, 3);
plot(tVec, jointAngles(:, 3), 'LineWidth', 1.3);
grid on;
xlabel('Time (s)');
ylabel('Joint 3 angle (rad)');
title('Joint 3');

subplot(2, 2, 4);
plot(tVec, beamTipY, 'LineWidth', 1.4);
grid on;
xlabel('Time (s)');
ylabel('Beam tip y displacement (m)');
title('Flexible-base tip');

figurePath = fullfile(dataDir, 'joint_angles_and_beam_tip.png');
saveas(gcf, figurePath);
fprintf('Flexible-base manipulator simulation completed.\n');
fprintf('  max |joint angles| = [%s] rad\n', ...
    sprintf('%.6g ', max(abs(jointAngles), [], 1)));
fprintf('  max |beam tip y| = %.6e m\n', max(abs(beamTipY)));
fprintf('  data: %s\n', outputCsv);
fprintf('  figure: %s\n', figurePath);

function dx = flexible_base_rhs(~, x, model, NB, nEta, jointTorque)
q = x(1:NB);
qd = x(NB + (1:NB));
eta = x(2 * NB + (1:nEta));
eta_d = x(2 * NB + nEta + (1:nEta));

tau = [0; 0; jointTorque; jointTorque; jointTorque];
[qdd, eta_dd] = FD_flex(model, q, qd, tau, eta, eta_d);
dx = [qd; qdd; eta_d; eta_dd{1}(:)];
end
