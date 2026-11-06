clear; close all; clc;
clear classes;

projectDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(projectDir, 'src', 'classes'), ...
        fullfile(projectDir, 'src', 'dynamics'), fullfile(projectDir, 'src', 'kinematics'), ...
        fullfile(projectDir, 'src', 'utils'), ...
        fullfile(projectDir, 'src', 'tests', 'dynamics', 'support'));

assert(~isempty(which('casadi.MX')), '请先将 CasADi 添加到 MATLAB 路径。');

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
nArm = 7;
armLengths = [0.28; 0.30; 0.26; 0.22; 0.18; 0.14; 0.10];
jointTypes = {'Rz'; 'Ry'; 'Ry'; 'Rx'; 'Ry'; 'Rx'; 'Ry'};

armBodies = cell(nArm, 1);
for i = 1:nArm
    p.mass = 0.8 + 0.1*i;
    p.inertia = [0.02+0.003*i; 0.03+0.004*i; 0.04+0.005*i];
    p.com = [0.5*armLengths(i); 0; 0];
    armBodies{i} = RigidBody(p);
end

NB = 8;
model.NB = NB;
model.parent = [0; 1; 2; 3; 4; 5; 6; 7];
model.jtype = [jointTypes; {'Rz'}];
model.Xtree = cell(NB, 1);
model.Xtree{1} = eye(6);

for i = 2:nArm
    model.Xtree{i} = xlt([armLengths(i-1), 0, 0]);
end
model.Xtree{8} = xlt([armLengths(end), 0, 0]);

model.rigidbody = [armBodies; {[]}];
model.flexbody = [repmat({[]}, nArm, 1); {rod}];
model.link_type = [repmat({'rigid'}, nArm, 1); {'flexible'}];
model.fixed_joints = 8;
model.gravity = zeros(3, 1);

model.Phi = repmat({zeros(6, 0)}, NB, 1);
model.Phi{8} = zeros(6, nEta);

model.modal_idx = repmat({[]}, NB, 1);
model.modal_idx{8} = 1:nEta;

%% Constant external forces
f_ext = cell(NB, 1);
f_modal = cell(NB, 1);

%% Build automatic-differentiation function using symbolic variables
% Keep symbolic inputs separate from the numeric state used by ID_flex_FO.
engine = 'MX';

q_sym   = casadi.MX.sym('q',   NB,   1);
qd_sym  = casadi.MX.sym('qd',  NB,   1);
qdd_sym = casadi.MX.sym('qdd', NB,   1);
eta_sym = casadi.MX.sym('eta', nEta, 1);
ed_sym  = casadi.MX.sym('ed',  nEta, 1);
edd_sym = casadi.MX.sym('edd', nEta, 1);

[tau_sym, tf_sym] = ID_flex( ...
    model, q_sym, qd_sym, qdd_sym, ...
    eta_sym, ed_sym, edd_sym, f_ext, f_modal);

% Rows: [all NB joint torques; all modal forces].
force_sym = [tau_sym; vertcat(tf_sym{:})];

% Columns: [q; eta] and [qd; ed].
tau_q_sym  = jacobian(force_sym, [q_sym; eta_sym]);
tau_qd_sym = jacobian(force_sym, [qd_sym; ed_sym]);

AD = casadi.Function(['flex_state_AD_', engine], ...
    {q_sym, qd_sym, qdd_sym, eta_sym, ed_sym, edd_sym}, ...
    {tau_q_sym, tau_qd_sym});

%% Evaluate analytic and AD Jacobians at the same numeric state
rng(1);

q   = 0.2*randn(NB, 1);
qd  = 0.3*randn(NB, 1);
qdd = 0.4*randn(NB, 1);
eta = 1e-5*randn(nEta, 1);
ed  = 1e-4*randn(nEta, 1);
edd = 1e-3*randn(nEta, 1);

% The eighth joint is fixed.
q(model.fixed_joints) = 0;
qd(model.fixed_joints) = 0;
qdd(model.fixed_joints) = 0;

[tau_q, tau_qd] = ID_flex_FO( ...
    model, q, qd, qdd, eta, ed, edd, f_ext, f_modal);

[tau_q_AD, tau_qd_AD] = AD(q, qd, qdd, eta, ed, edd);
tau_q_AD = full(tau_q_AD);
tau_qd_AD = full(tau_qd_AD);

%% Compare results
error_q = tau_q - tau_q_AD;
error_qd = tau_qd - tau_qd_AD;

maxError_q = max(abs(error_q(:)));
maxError_qd = max(abs(error_qd(:)));

relativeError_q = norm(error_q, 'fro') / max(1, norm(tau_q_AD, 'fro'));
relativeError_qd = norm(error_qd, 'fro') / max(1, norm(tau_qd_AD, 'fro'));

fprintf('Jacobian size: %d x %d\n', size(tau_q, 1), size(tau_q, 2));
fprintf('Position Jacobian: max abs = %.6e, relative = %.6e\n', ...
    maxError_q, relativeError_q);
fprintf('Velocity Jacobian: max abs = %.6e, relative = %.6e\n', ...
    maxError_qd, relativeError_qd);

assert(relativeError_q < 1e-9, '位置雅可比与 AD 结果不一致。');
assert(relativeError_qd < 1e-9, '速度雅可比与 AD 结果不一致。');
