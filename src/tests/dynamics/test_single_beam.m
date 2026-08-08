%% Single-link Flexible vs. Rigid Trajectory Comparison
clear; close all; clc;

% 添加路径
thisFile = mfilename('fullpath');
testDir = fileparts(thisFile);
srcDir = fileparts(fileparts(testDir));
projectDir = fileparts(srcDir);
refDir = fullfile(projectDir, 'refs', 'spatial_v2.1');
addpath(genpath(refDir), fullfile(srcDir, 'classes'), ...
        fullfile(srcDir, 'dynamics'), fullfile(srcDir, 'utils'));

%% 1. Beam Geometry and Physical Parameters
Lbeam  = 2.0;
W_beam = 0.1;
T_beam = 0.05;
rho    = 2700.0;
E      = 7.0e7;
nu     = 0.33;
nElements = 2;

beam.length = Lbeam;
beam.width = W_beam;
beam.thickness = T_beam;
beam.density = rho;
beam.youngs_modulus = E;
beam.poisson_ratio = nu;
beam.elements_num = nElements;
beam.fixed_root = true;

m_total = Lbeam * W_beam * T_beam * rho; % 27 kg
fprintf('Beam total mass: %.2f kg\n', m_total);

%% 2. Flexible Model Construction (FD_flex)
femBeam = FEM_Beam(beam);
[~, K_beam] = femBeam.get_matrices_3d();
[feShapeMatrix, stiffnessDofs] = femBeam.build_direct_fe_shape_matrix();
nFlexibleDofs = size(feShapeMatrix, 2);

Kff = full(K_beam(stiffnessDofs, stiffnessDofs));
Kff = 0.5 * (Kff + Kff.');

% FFRF 离散积分点
nx = 4; ny = 3; nz = 3;
[u0, S, weights] = femBeam.build_ffrf_data(feShapeMatrix, nx, ny, nz);

% 构建 FlexBody 对象
Dff = zeros(nFlexibleDofs);
params.density = rho;
params.u0      = u0;
params.S       = S;
params.weights = weights;
params.Kff     = Kff;
params.Dff     = Dff;
params.q       = [1; 0; 0; 0; 0; 0; 0; zeros(nFlexibleDofs, 1)];
params.qdot    = zeros(7 + nFlexibleDofs, 1);
rod = FlexBody(params);

% 构建柔性关节模型
NB = 1;
flexModel.NB = NB;
flexModel.parent = 0;
flexModel.jtype = {'Rz'};
flexModel.Xtree = {eye(6)};
flexModel.flexbody = {rod};
flexModel.link_type = {'flexible'};

% 提取末端 6D 空间形貌矩阵 [转角; 线位移]
tipNodeId = nElements + 1;
tipDofs = femBeam.get_node_dofs(tipNodeId);
tipShape6D = feShapeMatrix(tipDofs, :);
Phi_tip = [tipShape6D(4:6, :); tipShape6D(1:3, :)];
flexModel.Phi = {Phi_tip};
flexModel.gravity = zeros(6, 1);

%% 3. Equivalent Rigid Model Construction (spatial_v2 / FDab)
% 计算质心惯量 (Centroidal Inertia Tensor)
I_xx = (1/12) * m_total * (W_beam^2 + T_beam^2);
I_yy = (1/12) * m_total * (Lbeam^2 + T_beam^2);
I_zz = (1/12) * m_total * (Lbeam^2 + W_beam^2);
I_com = diag([I_xx, I_yy, I_zz]);

rigidParams.mass = m_total;
rigidParams.com = [Lbeam/2; 0; 0]; % 质心位于梁中心
rigidParams.inertia = I_com;
rigidBodyObj = RigidBody(rigidParams);

rigidModel.NB = NB;
rigidModel.parent = 0;
rigidModel.jtype = {'Rz'};
rigidModel.Xtree = {eye(6)};
rigidModel.gravity = zeros(6, 1);
rigidModel.I = {rigidBodyObj.get_spatial_inertia()};

%% 4. Dynamic Time Integration (Semi-implicit Euler)
dt = 1e-4;         % 步长 (柔性体高频振荡需较小步长)
T = 20.0;           % 仿真总时长
N = round(T / dt);

% 柔性体初始状态
qFlex = 0;
qdFlex = 0;
eta = zeros(nFlexibleDofs, 1);
eta_d = zeros(nFlexibleDofs, 1);

% 刚体初始状态
qRigid = 0;
qdRigid = 0;

% 历史记录开辟空间
qHistFlex  = zeros(1, N + 1);
qHistRigid = zeros(1, N + 1);
qdHistFlex = zeros(1, N + 1);
qdHistRigid= zeros(1, N + 1);
tVec       = (0:N) * dt;

for k = 1:N
    t = (k - 1) * dt;
    tau = 1.0; % 恒定驱动力矩 (也可改为正弦: 0.1*sin(2*pi*t))
    
    % 柔性体前向动力学
    [qddFlex, eta_dd_cell] = FD_flex(flexModel, qFlex, qdFlex, tau, eta, eta_d);
    eta_dd = eta_dd_cell{1};
    
    % 刚体前向动力学 (spatial_v2)
    qddRigid = FDab(rigidModel, qRigid, qdRigid, tau);
    
    % 半隐式欧拉积分 (Semi-implicit Euler)
    qdFlex  = qdFlex + dt * qddFlex;
    qFlex   = qFlex + dt * qdFlex;
    eta_d   = eta_d + dt * eta_dd;
    eta     = eta + dt * eta_d;
    
    qdRigid = qdRigid + dt * qddRigid;
    qRigid  = qRigid + dt * qdRigid;
    
    qHistFlex(k + 1)  = qFlex;
    qHistRigid(k + 1) = qRigid;
    qdHistFlex(k + 1) = qdFlex;
    qdHistRigid(k + 1)= qdRigid;
end

%% 5. Visualization and Analysis
figure('Color', 'w', 'Name', 'Flexible vs Rigid Link Comparison');

subplot(2, 1, 1);
plot(tVec, qHistFlex, 'r-', 'LineWidth', 1.2); hold on;
plot(tVec, qHistRigid, 'k--', 'LineWidth', 1.4);
grid on; ylabel('Joint Angle q (rad)');
legend('Flexible Beam (FD\_flex)', 'Rigid Beam (spatial\_v2)', 'Location', 'northwest');
title('Joint Trajectory Comparison');

subplot(2, 1, 2);
plot(tVec, qHistFlex - qHistRigid, 'b-', 'LineWidth', 1.0);
grid on; xlabel('Time (s)'); ylabel('Error: q_{flex} - q_{rigid} (rad)');
title('Position Deviation due to Structural Flexibility');

fprintf('Simulation finished. Max deviation = %.5e rad\n', max(abs(qHistFlex - qHistRigid)));