%% Seven-DOF forward-dynamics trajectory comparison
% Compares FD_flex with the spatial_v2.1 FDab implementation under the
% same time-varying joint torques.

thisFile = mfilename('fullpath');
testDir = fileparts(thisFile);
srcDir = fileparts(fileparts(testDir));
projectDir = fileparts(srcDir);
refDir = fullfile(projectDir, 'refs', 'spatial_v2.1');
addpath(genpath(refDir), fullfile(srcDir, 'classes'), ...
        fullfile(srcDir, 'dynamics'), fullfile(srcDir, 'kinematics'), fullfile(srcDir, 'utils'));

n = 7;
body = cell(n, 1);
for i = 1:n
    p.mass = 0.7 + 0.08 * i;
    p.inertia = [0.03 + 0.005*i; 0.05 + 0.006*i; 0.08 + 0.008*i];
    p.com = [0; 0; 0];
    body{i} = RigidBody(p);
end

model.NB = n;
model.parent = [0; 1; 2; 3; 4; 5; 6];
model.jtype = {'Rz'; 'Ry'; 'Ry'; 'Rx'; 'Ry'; 'Rx'; 'Ry'};
model.Xtree = cell(n, 1);
for i = 1:n
    model.Xtree{i} = xlt([0.25 + 0.04*i, 0, 0]);
end
model.rigidbody = body;
model.link_type = repmat({'rigid'}, n, 1);
model.Phi = repmat({zeros(6, 0)}, n, 1);
model.gravity = [0; 0; -9.81];

refModel.NB = n;
refModel.parent = model.parent;
refModel.jtype = model.jtype;
refModel.Xtree = model.Xtree;
refModel.gravity = model.gravity;
refModel.I = cellfun(@(b) b.get_spatial_inertia(), body, 'UniformOutput', false);

qFlex = [0.15; -0.25; 0.30; -0.20; 0.10; 0.35; -0.18];
qdFlex = zeros(n, 1);
qRef = qFlex;
qdRef = qdFlex;

dt = 1e-3;
T = 5.0;
N = round(T / dt);
qHistoryFlex = zeros(n, N + 1);
qHistoryRef = zeros(n, N + 1);
qHistoryFlex(:, 1) = qFlex;
qHistoryRef(:, 1) = qRef;

for k = 1:N
    t = (k - 1) * dt;
    tau = [0.30*sin(2*t); 0.20*cos(1.5*t); 0.15*sin(2.5*t); ...
           0.12*cos(1.2*t); 0.10*sin(1.8*t); 0.08*cos(2.2*t); ...
           0.05*sin(1.1*t)];

    qddFlex = FD_flex(model, qFlex, qdFlex, tau, [], []);
    qddRef = FDab(refModel, qRef, qdRef, tau);

    % Identical semi-implicit Euler integration for both models.
    qdFlex = qdFlex + dt * qddFlex;
    qFlex = qFlex + dt * qdFlex;
    qdRef = qdRef + dt * qddRef;
    qRef = qRef + dt * qdRef;
    qHistoryFlex(:, k + 1) = qFlex;
    qHistoryRef(:, k + 1) = qRef;
end

err = max(abs(qHistoryFlex - qHistoryRef), [], 'all');
fprintf('FD_flex 7-DOF trajectory comparison:\n');
fprintf('  max joint-position error = %.3e rad\n', err);
assert(err < 1e-9, 'FD_flex and spatial_v2 trajectories are inconsistent.');
fprintf('Trajectory comparison passed.\n');

figure('Name', 'FD_flex vs spatial_v2 - 7DOF');
plot((0:N)*dt, qHistoryFlex - qHistoryRef, 'LineWidth', 1.0);
xlabel('Time (s)'); ylabel('q_{FD\_flex} - q_{spatial\_v2} (rad)');
grid on; legend('q_1','q_2','q_3','q_4','q_5','q_6','q_7', 'Location', 'best');
