%% ID_flex test: FlexBody rigid link connected to a FlexBody flexible beam
% Run this file directly from the project root or from any current folder.
% No model.I or manually supplied inertia blocks are used.

thisFile = mfilename('fullpath');
testDir = fileparts(thisFile);
srcDir = fileparts(fileparts(testDir));
supportDir = fullfile(testDir, 'support');
addpath(fullfile(srcDir, 'classes'), fullfile(srcDir, 'dynamics'), fullfile(srcDir, 'kinematics'), ...
        fullfile(srcDir, 'utils'), supportDir);

%% Build a rigid link with directly specified inertia
rigidParams.mass = 2;
rigidParams.inertia = [0.12; 0.18; 0.25];
rigidParams.com = [0; 0; 0];
rigidBody = RigidBody(rigidParams);

%% Build a flexible beam (one translational modal coordinate)
beamParams.u0 = [1; 0; 0];
beamParams.S = reshape([0; 0; 1], 3, 1, 1);
beamParams.massWeights = 5;
beamParams.Kff = 7;
beamParams.Dff = 2;
beamParams.q = [1; 0; 0; 0; 0; 0; 0; 0.1];
beamParams.qdot = [zeros(7, 1); 0.2];
beamBody = FlexBody(beamParams);

%% Serial topology: rigid Rz link -> flexible Rz beam
model.NB = 2;
model.parent = [0; 1];
model.jtype = {'Rz'; 'Rz'};
model.Xtree = {eye(6); eye(6)};
model.rigidbody = {rigidBody; []};
model.flexbody = {[]; beamBody};
model.link_type = {'rigid'; 'flexible'};
model.Phi = {zeros(6, 0); [0; 0; 0; 0; 0; 1]};
model.gravity = zeros(6, 1);

q = [0; 0];
qd = [0; 0];
qdd = [1.2; 2.0];
eta = 0.1;
eta_d = 0.2;

%% Run inverse dynamics
[tau, eta_dd] = ID_flex(model, q, qd, qdd, eta, eta_d);

assert(isequal(size(tau), [2, 1]), ...
    'ID_flex must return one torque per active joint.');
assert(iscell(eta_dd) && numel(eta_dd) == 2 && isempty(eta_dd{1}) ...
       && isequal(size(eta_dd{2}), [1, 1]), ...
    'The flexible beam modal acceleration was not returned correctly.');
assert(all(isfinite(tau)) && all(isfinite(eta_dd{2})), ...
    'ID_flex returned NaN or Inf.');

fprintf('FlexBody-only ID_flex test passed.\n');
fprintf('  tau    = [%+.12g, %+.12g]^T\n', tau(1), tau(2));
fprintf('  eta_dd = %+.12g\n', eta_dd{2});
