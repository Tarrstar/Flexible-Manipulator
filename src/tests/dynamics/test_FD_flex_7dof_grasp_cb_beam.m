%% 7-DOF arm grasping a Craig-Bampton flexible rod
clear; close all; clc;
clear classes;
thisFile = mfilename('fullpath');
testDir = fileparts(thisFile);
srcDir = fileparts(fileparts(testDir));
projectDir = fileparts(srcDir);
dataDir = fullfile(projectDir, 'data', 'arm_grasp_cb_beam');
if ~exist(dataDir, 'dir'), mkdir(dataDir); end
addpath(fullfile(srcDir, 'classes'), fullfile(srcDir, 'dynamics'), ...
        fullfile(srcDir, 'utils'), fullfile(testDir, 'support'));

%% CB-reduced flexible rod
beam.length = 1.0; 
beam.width = 0.04;
beam.thickness = 0.025;
beam.density = 2700; 
beam.youngs_modulus = 7e7; 
beam.poisson_ratio = 0.33;
beam.elements_num = 4;
beam.fixed_root = true; 
beam.nInteriorModes = 18;
beam.boundaryNodeIds = [1, beam.elements_num + 1];
femBeam = FEM_Beam_CB_Reduction(beam);
[u0, Sbeam, weights, cbShape, Kff] = ...
    femBeam.build_ffrf_data(beam.nInteriorModes, 4, 3, 3);
nEta = size(cbShape, 2);
rodParams.density = beam.density;
rodParams.u0 = u0; 
rodParams.S = Sbeam; 
rodParams.weights = weights;
rodParams.Kff = Kff;
rodParams.Dff = zeros(nEta);
rodParams.q = [1; 0; 0; 0; zeros(3+nEta, 1)];
rodParams.qdot = zeros(7+nEta, 1);
rod = FlexBody(rodParams);

%% Seven-DOF arm and fixed grasp connection
nArm = 7; armLengths = [0.28; 0.30; 0.26; 0.22; 0.18; 0.14; 0.10];
jointTypes = {'Rz'; 'Ry'; 'Ry'; 'Rx'; 'Ry'; 'Rx'; 'Ry'};
armBodies = cell(nArm, 1);
for i = 1:nArm
    p.mass = 0.8 + 0.1*i;
    p.inertia = [0.02+0.003*i; 0.03+0.004*i; 0.04+0.005*i];
    p.com = [0.5 * armLengths(i); 0; 0];
    armBodies{i} = RigidBody(p);
end

% Link 8 is virtual and fixed; only seven arm torques are prescribed.
NB = 8; 
model.NB = NB;
model.parent = [0; 1; 2; 3; 4; 5; 6; 7];
model.jtype = [jointTypes; {'Rz'}];
model.Xtree = cell(NB, 1);
model.Xtree{1} = eye(6);
for i = 2:nArm
    model.Xtree{i} = xlt([armLengths(i - 1), 0, 0]);
end
model.Xtree{8} = xlt([armLengths(end), 0, 0]);
model.rigidbody = [armBodies; {[]}];
model.flexbody = [{[]}; {[]}; {[]}; {[]}; {[]}; {[]}; {[]}; {rod}];
model.link_type = [repmat({'rigid'}, nArm, 1); {'flexible'}];
model.fixed_joints = 8;
tipShape = femBeam.evaluate_displacement(cbShape, femBeam.length, 0, 0);
model.Phi = repmat({zeros(6, 0)}, NB, 1);
model.Phi{8} = zeros(6, nEta);
% model.gravity = [0; 0; -9.81];

%% Torque-driven simulation
T = 1.0;
flexLink = find(strcmpi(model.link_type, 'flexible'), 1, 'first');
if isempty(flexLink)
    error('The test model must contain one flexible link.');
end

% State vector: [q; qd; eta; eta_d].  ode45 uses adaptive time steps.
x0 = zeros(2*NB + 2*nEta, 1);
odeOptions = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[tVec, xHist] = ode45(@(t, x) grasp_cb_rhs(t, x, model, NB, nEta, flexLink), ...
                      [0, T], x0, odeOptions);

qHist = xHist(:, 1:NB);
etaHist = xHist(:, 2*NB + (1:nEta));
tipY = zeros(numel(tVec), 1);
graspY = zeros(numel(tVec), 1);
for k = 1:numel(tVec)
    [tipY(k), graspY(k)] = grasped_tip_y(model, qHist(k, :).', ...
        etaHist(k, :).', femBeam.length, tipShape);
end

matlabCsv = fullfile(dataDir, 'matlab_cb_tip_y.csv');
writematrix_with_header(matlabCsv, [tVec(:), tipY, graspY], ...
    'time,matlab_tip_y,matlab_grasp_y');

%% Compare with Exudyn 3D beam endpoint-y export
exudynCsv = fullfile(dataDir, 'exudyn_beam_tip_y.csv');
exuY = [];
exuGraspY = [];
if exist(exudynCsv, 'file')
    exu = readmatrix(exudynCsv); exu = exu(all(isfinite(exu), 2), :);
    [~, uniqueIdx] = unique(exu(:, 1), 'stable');
    exu = exu(uniqueIdx, :);
    if size(exu, 2) >= 3
        exuGraspY = interp1(exu(:,1), exu(:,3), tVec(:), ...
            'linear', 'extrap');
        graspErr = graspY - exuGraspY;
    else
        exuGraspY = [];
        graspErr = [];
    end
    if size(exu, 2) >= 8, exu = exu(:, [1, 7]);
    elseif size(exu, 2) >= 2, exu = exu(:, 1:2);
    else, error('Exudyn CSV needs time and tip-y columns.'); end
    exuY = interp1(exu(:,1), exu(:,2), tVec(:), 'linear', 'extrap');
    err = tipY - exuY;
    fprintf('CB-FlexBody vs Exudyn 3D beam endpoint-y comparison:\n');
    fprintf('  RMS error = %.6e m\n', sqrt(mean(err.^2)));
    fprintf('  max error = %.6e m\n', max(abs(err)));
    if ~isempty(graspErr)
        fprintf('  grasp RMS error = %.6e m\n', ...
            sqrt(mean(graspErr.^2)));
        fprintf('  grasp max error = %.6e m\n', max(abs(graspErr)));
    end
else
    warning('Exudyn data not found: %s', exudynCsv); exuY = [];
end

%% Compare with the Exudyn equal-mass, equal-inertia rigid beam
exudynRigidCsv = fullfile(dataDir, 'exudyn_rigid_beam_tip_y.csv');
exuRigidY = [];
exuRigidGraspY = [];
if exist(exudynRigidCsv, 'file')
    exuRigid = readmatrix(exudynRigidCsv);
    exuRigid = exuRigid(all(isfinite(exuRigid), 2), :);
    if size(exuRigid, 2) < 2
        error('Exudyn rigid-beam CSV needs time and tip-y columns.');
    end
    [~, uniqueIdx] = unique(exuRigid(:, 1), 'stable');
    exuRigid = exuRigid(uniqueIdx, :);
    exuRigidY = interp1(exuRigid(:, 1), exuRigid(:, 2), ...
        tVec(:), 'linear', 'extrap');
    rigidTipErr = tipY - exuRigidY;
    fprintf('CB-FlexBody vs Exudyn equivalent rigid beam comparison:\n');
    fprintf('  tip RMS difference = %.6e m\n', ...
        sqrt(mean(rigidTipErr.^2)));
    fprintf('  tip max difference = %.6e m\n', max(abs(rigidTipErr)));
    if size(exuRigid, 2) >= 3
        exuRigidGraspY = interp1(exuRigid(:, 1), exuRigid(:, 3), ...
            tVec(:), 'linear', 'extrap');
        rigidGraspErr = graspY - exuRigidGraspY;
        fprintf('  grasp RMS difference = %.6e m\n', ...
            sqrt(mean(rigidGraspErr.^2)));
        fprintf('  grasp max difference = %.6e m\n', ...
            max(abs(rigidGraspErr)));
    end
else
    warning(['Run refs/exudyn_7dof_grasp_rigid_beam.py first: %s'], ...
        exudynRigidCsv);
end

figure('Color', 'w', 'Name', '7DOF grasped beam comparison');
subplot(2, 1, 1);
plot(tVec, graspY, 'LineWidth', 1.4); hold on;
graspLegend = {'MATLAB CB-FlexBody'};
if ~isempty(exuGraspY)
    plot(tVec, exuGraspY, '--', 'LineWidth', 1.2);
    graspLegend{end + 1} = 'Exudyn 3D beam';
end
if ~isempty(exuRigidGraspY)
    plot(tVec, exuRigidGraspY, ':', 'LineWidth', 1.5);
    graspLegend{end + 1} = 'Exudyn equivalent rigid beam';
end
grid on; ylabel('Grasp point y (m)');
legend(graspLegend, 'Location', 'best');

subplot(2, 1, 2);
plot(tVec, tipY, 'LineWidth', 1.4); hold on;
tipLegend = {'MATLAB CB-FlexBody'};
if ~isempty(exuY)
    plot(tVec, exuY, '--', 'LineWidth', 1.2);
    tipLegend{end + 1} = 'Exudyn 3D beam';
end
if ~isempty(exuRigidY)
    plot(tVec, exuRigidY, ':', 'LineWidth', 1.5);
    tipLegend{end + 1} = 'Exudyn equivalent rigid beam';
end
grid on; xlabel('Time (s)'); ylabel('Rod endpoint y (m)');
legend(tipLegend, 'Location', 'best');
saveas(gcf, fullfile(dataDir, 'cb_beam_tip_y_vs_exudyn.png'));
fprintf('MATLAB data saved to %s\n', matlabCsv);

function [y, graspY] = grasped_tip_y(model, q, eta, beamLength, tipShape)
X = eye(6);
for i = 1:model.NB
    [XJ, ~] = jcalc(model.jtype{i}, q(i));
    X = XJ * model.Xtree{i} * X;
end
tipLocal = [beamLength; 0; 0] + tipShape * eta;
tip = Xpt(inv(X), tipLocal);
y = tip(2);
grasp = Xpt(inv(X), zeros(3, 1));
graspY = grasp(2);
end

function dx = grasp_cb_rhs(t, x, model, NB, nEta, flexLink)
q = x(1:NB);
qd = x(NB + (1:NB));
eta = x(2*NB + (1:nEta));
eta_d = x(2*NB + nEta + (1:nEta));

tau = [0.35*sin(2*pi*0.8*t); 0.25*cos(2*pi*0.6*t); ...
       0.20*sin(2*pi*0.9*t); 0.16*cos(2*pi*0.7*t); ...
       0.12*sin(2*pi*0.5*t); 0.10*cos(2*pi*0.8*t); ...
       0.08*sin(2*pi*0.4*t); 0];
[qdd, eta_dd] = FD_flex(model, q, qd, tau, eta, eta_d);

dx = [qd; qdd; eta_d; eta_dd{flexLink}(:)];
end
