%% Planar rigid-flexible double pendulum: MATLAB FD_flex vs spatial_v2 Rigid Reference
clear; close all; clc;
clear classes;
thisFile = mfilename('fullpath');
testDir = fileparts(thisFile);
srcDir = fileparts(fileparts(testDir));
projectDir = fileparts(srcDir);
dataDir = fullfile(projectDir, 'data', 'planar_rigid_flex');
if ~exist(dataDir, 'dir'), mkdir(dataDir); end

addpath(fullfile(srcDir, 'classes'), fullfile(srcDir, 'dynamics'), ...
        fullfile(srcDir, 'utils'), fullfile(testDir, 'support'), ...
        fullfile(projectDir, 'refs', 'spatial_v2.1', 'spatial'));

%% 1. Geometry and Flexible Beam (Craig-Bampton Reduction)
L1 = 0.50;
L2 = 0.50;

beam.length = L2;
beam.width = 0.04;
beam.thickness = 0.025;
beam.density = 2700;
beam.youngs_modulus = 7e7;
beam.poisson_ratio = 0.33;
beam.elements_num = 4;
beam.fixed_root = true;
beam.boundaryNodeIds = [1, beam.elements_num + 1];
beam.nInteriorModes = 12;

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

%% 2. Flexible Model Setup (FD_flex)
NB = 2;
flexModel.NB = NB;
flexModel.parent = [0; 1];
flexModel.jtype = {'Rz'; 'Rz'};
flexModel.Xtree = {eye(6), xlt([L1, 0, 0])};
flexModel.rigidbody = cell(NB, 1);
flexModel.flexbody = cell(NB, 1);

rigidParams.mass = 1.0;
rigidParams.inertia = diag([0.01, 0.01, 1.0 * L1^2 / 12]);
rigidParams.com = [L1/2; 0; 0];
flexModel.rigidbody{1} = RigidBody(rigidParams);
flexModel.flexbody{2} = rod;
flexModel.link_type = {'rigid'; 'flexible'};
flexModel.Phi = {zeros(6, 0), zeros(6, nEta)};
flexModel.gravity = zeros(6, 1);

%% 3. Pure Rigid Model Setup (spatial_v2 / FDab)
rigidModel.NB = NB;
rigidModel.parent = [0; 1];
rigidModel.jtype = {'Rz'; 'Rz'};
rigidModel.Xtree = {eye(6), xlt([L1, 0, 0])};
rigidModel.gravity = zeros(6, 1);

% Link 1 惯量 (Spatial inertia)
m1 = 1.0;
com1 = [L1/2; 0; 0];
I1_com = diag([0.01, 0.01, 1.0 * L1^2 / 12]);

% Link 2 惯量 (将柔性梁退化为刚体的转动惯量)
m2 = beam.density * (beam.width * beam.thickness) * L2; % 1.35 kg
com2 = [L2/2; 0; 0];
I2_x = (1/12) * m2 * (beam.width^2 + beam.thickness^2);
I2_y = (1/12) * m2 * (L2^2 + beam.thickness^2);
I2_z = (1/12) * m2 * (L2^2 + beam.width^2);
I2_com = diag([I2_x, I2_y, I2_z]);

rigidModel.I = {mcI(m1, com1, I1_com), mcI(m2, com2, I2_com)};

%% 4. Time Domain Simulations (ode45)
T = 1.0;
odeOptions = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

% --- 4.1 刚柔混合模型仿真 ---
x0_flex = zeros(2*NB + 2*nEta, 1);
[tVec, xHistFlex] = ode45(@(t, x) planar_rhs(t, x, flexModel, NB, nEta), ...
                          [0, T], x0_flex, odeOptions);

qHistFlex = xHistFlex(:, 1:NB);
etaHistFlex = xHistFlex(:, 2*NB + (1:nEta));
tipShape = femBeam.evaluate_displacement(cbShape, L2, 0, 0);
tipY_flex = zeros(numel(tVec), 1);
for k = 1:numel(tVec)
    tipY_flex(k) = planar_tip_y(flexModel, qHistFlex(k, :).', ...
        etaHistFlex(k, :).', L2, tipShape);
end

% --- 4.2 MATLAB spatial_v2 纯刚体模型仿真 ---
x0_rigid = zeros(2*NB, 1);
[tRigid, xHistRigid] = ode45(@(t, x) planar_rigid_rhs(t, x, rigidModel, NB), ...
                             tVec, x0_rigid, odeOptions);

qHistRigid = xHistRigid(:, 1:NB);
tipY_rigid = zeros(numel(tRigid), 1);
for k = 1:numel(tRigid)
    tipY_rigid(k) = planar_rigid_tip_y(qHistRigid(k, :).', L1, L2);
end

% 保存 MATLAB 采样数据
matlabCsv = fullfile(dataDir, 'matlab_endpoint_y.csv');
writematrix_with_header(matlabCsv, [tVec(:), tipY_flex, tipY_rigid], ...
                        'time,matlab_flex_tip_y,matlab_rigid_tip_y');

%% 5. Read Exudyn ANCF Data (if available)
exudynCsv = fullfile(dataDir, 'exudyn_endpoint_y.csv');
exuY = [];
if exist(exudynCsv, 'file')
    exu = readmatrix(exudynCsv);
    exu = exu(all(isfinite(exu), 2), :);
    [~, uIdx] = unique(exu(:,1));
    exu = exu(uIdx, :);
    exuY = interp1(exu(:,1), exu(:,2), tVec(:), 'linear', 'extrap');
    
    err = tipY_flex - exuY;
    fprintf('Planar rigid-flexible (FD_flex) vs Exudyn ANCF:\n');
    fprintf('  RMS error = %.6e m\n', sqrt(mean(err.^2)));
    fprintf('  max error = %.6e m\n', max(abs(err)));
else
    warning('Exudyn ANCF reference not found: %s', exudynCsv);
end

%% 6. Visualization (对比图)
figure('Color', 'w', 'Name', 'Planar Double Pendulum Comparison');
plot(tVec, tipY_flex, 'r-', 'LineWidth', 1.5); hold on;
plot(tRigid, tipY_rigid, 'k-.', 'LineWidth', 1.4);

legendEntries = {'MATLAB Rigid-Flexible (FD\_flex)', ...
                 'MATLAB Pure Rigid (spatial\_v2)'};

if ~isempty(exuY)
    plot(tVec, exuY, 'b--', 'LineWidth', 1.2);
    legendEntries{end+1} = 'Exudyn ANCF (Flexible)';
end

grid on; xlabel('Time (s)'); ylabel('Endpoint Y-displacement (m)');
legend(legendEntries, 'Location', 'best');
title('Planar Double Pendulum Tip Trajectory Comparison');

saveas(gcf, fullfile(dataDir, 'planar_rigid_flex_endpoint_y.png'));
fprintf('Simulation complete. Figure saved.\n');

%% Helper Functions

% 刚柔混合微分方程
function dx = planar_rhs(~, x, model, NB, nEta)
q = x(1:NB);
qd = x(NB + (1:NB));
eta = x(2*NB + (1:nEta));
eta_d = x(2*NB + nEta + (1:nEta));

tau = [0.1; 0.0];
[qdd, eta_dd] = FD_flex(model, q, qd, tau, eta, eta_d);

dx = [qd; qdd; eta_d; eta_dd{2}(:)];
end

% spatial_v2 纯刚体微分方程
function dx = planar_rigid_rhs(~, x, rigidModel, NB)
q = x(1:NB);
qd = x(NB + (1:NB));

tau = [0.1; 0.0];
qdd = FDab(rigidModel, q, qd, tau);

dx = [qd; qdd];
end

% 刚柔混合模型末端 Y 坐标计算
function y = planar_tip_y(model, q, eta, beamLength, tipShape)
X = eye(6);
for i = 1:model.NB
    [XJ, ~] = jcalc(model.jtype{i}, q(i));
    if strcmpi(model.link_type{i}, 'flexible')
        [R, ~, delta] = get_deformation_transform(model.Phi{i}, eta);
        Xdef = [R, zeros(3); -R*get_skew_sym_matrix(delta(4:6)), R];
    else
        Xdef = eye(6);
    end
    X = Xdef * XJ * model.Xtree{i} * X;
end
tipLocal = [beamLength; 0; 0] + tipShape * eta;
tip = Xpt(inv(X), tipLocal);
y = tip(2);
end

% 纯刚体模型末端 Y 坐标正运动学计算
function y = planar_rigid_tip_y(q, L1, L2)
theta1 = q(1);
theta2 = q(1) + q(2);
y = L1 * sin(theta1) + L2 * sin(theta2);
end
