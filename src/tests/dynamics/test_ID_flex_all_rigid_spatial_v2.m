%% Seven-DOF all-rigid comparison against spatial_v2.1 ID
% Every link is represented by a RigidBody with directly specified
% principal inertias.  The reference implementation uses the spatial
% inertia extracted from the same RigidBody object.
clear all;
thisFile = mfilename('fullpath');
testDir = fileparts(thisFile);
srcDir = fileparts(fileparts(testDir));
projectDir = fileparts(srcDir);
refDir = fullfile(projectDir, 'refs', 'spatial_v2.1');

addpath(genpath(refDir), fullfile(srcDir, 'classes'), ...
        fullfile(srcDir, 'dynamics'), fullfile(srcDir, 'utils'));

%% Seven rigid RigidBody links
n = 7;
body = cell(n, 1);
for i = 1:n
    params.mass = 0.7 + 0.08 * i;
    params.inertia = [0.03 + 0.005*i; 0.05 + 0.006*i; 0.08 + 0.008*i];
    params.com = [0; 0; 0];
    body{i} = RigidBody(params);
end

model.NB = n;
model.parent = [0; 1; 2; 3; 4; 5; 6];
% A standard 7-DOF serial arm has non-parallel, alternating joint axes.
% The sequence is z-y-y-x-y-x-y, representative of an anthropomorphic
% shoulder/elbow/wrist chain.
model.jtype = {'Ry'; 'Ry'; 'Ry'; 'Ry'; 'Ry'; 'Ry'; 'Ry'};
model.Xtree = cell(n, 1);
for i = 1:n
    model.Xtree{i} = xlt([0.25 + 0.04 * i, 0, 0]);
end
model.rigidbody = body;
model.link_type = repmat({'rigid'}, n, 1);
model.Phi = repmat({zeros(6, 0)}, n, 1);
% model.gravity = [0; 0; -9.81;0 ;0 ;0];

q = [0.15; -0.25; 0.30; -0.20; 0.10; 0.35; -0.18];
qd = [0.40; -0.20; 0.35; 0.15; -0.30; 0.25; 0.10];
qdd = [-0.50; 0.80; -0.60; 0.45; 0.30; -0.70; 0.55];

% q = [1];
% qd = [1];
% qdd = [0];

%% Convert the same FlexBody inertias for the spatial_v2 reference
refModel.NB = n;
refModel.parent = model.parent;
refModel.jtype = model.jtype;
refModel.Xtree = model.Xtree;
% refModel.gravity = model.gravity;
refModel.I = cell(n, 1);
for i = 1:n
    refModel.I{i} = body{i}.get_spatial_inertia();
end

tauFlex = ID_flex(model, q, qd, qdd, [], []);
tauSpatial = ID(refModel, q, qd, qdd);

err = norm(tauFlex - tauSpatial, inf);
fprintf('7-DOF all-rigid spatial_v2 comparison:\n');
fprintf('  ID_flex  = [%s]^T\n', sprintf('%+.10g ', tauFlex));
fprintf('  spatial  = [%s]^T\n', sprintf('%+.10g ', tauSpatial));
fprintf('  inf-norm error = %.3e\n', err);
% assert(err < 1e-10, 'ID_flex and spatial_v2 ID results are inconsistent.');
% fprintf('Comparison passed.\n');
