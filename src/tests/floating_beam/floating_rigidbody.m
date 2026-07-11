clear; clc; close all;

thisFile = mfilename('fullpath');
thisDir = fileparts(thisFile);
srcDir = fileparts(fileparts(thisDir));
projectDir = fileparts(srcDir);
dataDir = fullfile(projectDir, 'data', 'floating_beam');

if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

addpath(thisDir, fullfile(srcDir, 'classes'), fullfile(srcDir, 'utils'));

%% ================== 1. Rigid beam parameters ==================

Lbeam  = 5.0;
W_beam = 0.1;
T_beam = 0.05;

rho = 2700.0;

Fy = 10.0;
Mz = 0.0;

%% ================== 2. Rigid FlexBody data ==================

nx = 6;
ny = 3;
nz = 3;

[u0, S, weights] = build_rigid_beam_data( ...
    Lbeam, W_beam, T_beam, nx, ny, nz);

nFlexibleDofs = 0;

params.density = rho;
params.u0 = u0;
params.S = S;
params.weights = weights;
params.Kff = zeros(nFlexibleDofs);
params.Dff = zeros(nFlexibleDofs);

% q = [theta; R]
% R is the left/base material point, matching the beam reference origin.
params.q = [1; 0; 0; 0; ...
            0; 0; 0];
params.qdot = zeros(7, 1);

body = FlexBody(params);

nq = 7;

%% ================== 3. Time integration ==================

simulationTime = 2.0;
timeStep = 2e-3;

tVec = 0:timeStep:simulationTime;
x0 = [params.q; params.qdot];

odeOptions = odeset( ...
    'RelTol', 1e-8, ...
    'AbsTol', 1e-10, ...
    'MaxStep', timeStep);

rhs = @(t, x) rhs_rigidbody(t, x, body, Fy, Mz);
[~, xHist] = ode45(rhs, tVec, x0, odeOptions);

log_data = zeros(numel(tVec), 6);

for k = 1:numel(tVec)
    q = xHist(k, 1:nq).';
    qdot = xHist(k, nq+1:end).';

    [q, qdot] = normalize_quaternion_state(q, qdot);
    xHist(k, 1:nq) = q.';
    xHist(k, nq+1:end) = qdot.';

    [base_x, base_y, base_theta, tip_x, tip_y] = ...
        get_rigidbody_output(q, Lbeam);

    log_data(k, :) = [tVec(k), base_x, base_y, ...
                      base_theta, tip_x, tip_y];
end

%% ================== 4. Save and compare ==================

csvName = fullfile(dataDir, 'floating_rigidbody_flexbody.csv');
writematrix_with_header(csvName, log_data, ...
    'time,base_x,base_y,base_theta,tip_x,tip_y');

referenceCsvName = fullfile(dataDir, 'floating_rigid_beam_exudyn.csv');
comparisonFigName = fullfile(dataDir, ...
    'floating_rigidbody_flexbody_vs_exudyn.png');

reference_data = [];

if exist(referenceCsvName, 'file')
    reference_data = readmatrix(referenceCsvName);
    reference_data = reference_data(all(isfinite(reference_data), 2), :);

    if size(reference_data, 2) < 6
        warning('Reference CSV does not contain the expected 6 columns: %s', ...
            referenceCsvName);
        reference_data = [];
    end
else
    warning(['Rigid reference CSV was not found. Run ', ...
             'refs/floating_rigid_body.py first to create: %s'], ...
             referenceCsvName);
end

figure('Color', 'w', 'Position', [120, 100, 1050, 720]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(log_data(:, 1), log_data(:, 2), 'LineWidth', 1.2); hold on;
plot(log_data(:, 1), log_data(:, 3), 'LineWidth', 1.2);
baseLegend = {'FlexBody Base X', 'FlexBody Base Y'};

if ~isempty(reference_data)
    plot(reference_data(:, 1), reference_data(:, 2), '--', 'LineWidth', 1.2);
    plot(reference_data(:, 1), reference_data(:, 3), '--', 'LineWidth', 1.2);
    baseLegend = [baseLegend, {'Exudyn Base X', 'Exudyn Base Y'}];
end

grid on;
title('Rigid Body Base');
xlabel('Time / s');
ylabel('Base position');
legend(baseLegend{:}, 'Location', 'best');

nexttile;
plot(log_data(:, 1), log_data(:, 5), 'LineWidth', 1.2); hold on;
plot(log_data(:, 1), log_data(:, 6), 'LineWidth', 1.2);
tipLegend = {'FlexBody Tip X', 'FlexBody Tip Y'};

if ~isempty(reference_data)
    plot(reference_data(:, 1), reference_data(:, 5), '--', 'LineWidth', 1.2);
    plot(reference_data(:, 1), reference_data(:, 6), '--', 'LineWidth', 1.2);
    tipLegend = [tipLegend, {'Exudyn Tip X', 'Exudyn Tip Y'}];
end

grid on;
title('Rigid Body Tip');
xlabel('Time / s');
ylabel('Tip position');
legend(tipLegend{:}, 'Location', 'best');

saveas(gcf, comparisonFigName);

disp(['Rigid FlexBody data saved to ', csvName]);
disp(['Comparison figure saved to ', comparisonFigName]);

%% ========================================================================
% Local functions
% ========================================================================

function [u0, S, weights] = build_rigid_beam_data( ...
        L, W, T, nx, ny, nz)

    [xi_x, wi_x] = gauss_legendre_1d(nx);
    [xi_y, wi_y] = gauss_legendre_1d(ny);
    [xi_z, wi_z] = gauss_legendre_1d(nz);

    xNodes = L / 2 * (xi_x + 1);
    yNodes = W / 2 * xi_y;
    zNodes = T / 2 * xi_z;

    jac = L / 2 * W / 2 * T / 2;
    nGauss = nx * ny * nz;

    u0 = zeros(3, nGauss);
    S = zeros(3, 0, nGauss);
    weights = zeros(nGauss, 1);

    g = 0;

    for ix = 1:nx
        x = xNodes(ix);

        for iy = 1:ny
            y = yNodes(iy);

            for iz = 1:nz
                z = zNodes(iz);
                g = g + 1;

                u0(:, g) = [x; y; z];
                weights(g) = wi_x(ix) * wi_y(iy) * wi_z(iz) * jac;
            end
        end
    end
end

function dx = rhs_rigidbody(~, x, body, Fy, Mz)

    nq = numel(x) / 2;

    q = x(1:nq);
    qdot = x(nq+1:end);

    [q, qdot] = normalize_quaternion_state(q, qdot);

    M = body.get_mass_matrix(q);
    Qv = body.get_quadratic_velocity_vec(q, qdot);

    theta = q(1:4);
    thetaDot = qdot(1:4);

    [G, ~] = get_quaternion_G_matrices(theta);

    Qext = zeros(nq, 1);
    torque_global = [0; 0; Mz];
    Qext(1:4) = G.' * torque_global;
    Qext(5:7) = [0; Fy; 0];

    Phi_q = [2 * theta.', zeros(1, nq - 4)];
    gamma = -2 * (thetaDot.' * thetaDot);

    KKT = [M, -Phi_q.';
           Phi_q, 0];

    rhs = [Qext + Qv;
           gamma];

    sol = KKT \ rhs;
    qddot = sol(1:nq);

    dx = [qdot;
          qddot];
end

function [base_x, base_y, base_theta, tip_x, tip_y] = ...
        get_rigidbody_output(q, Lbeam)

    theta = q(1:4);
    R = q(5:7);

    A = get_rotation_matrix_quaternion(theta);

    base_x = R(1);
    base_y = R(2);
    base_theta = atan2(A(2, 1), A(1, 1));

    rTipGlobal = R + A * [Lbeam; 0; 0];

    tip_x = rTipGlobal(1);
    tip_y = rTipGlobal(2);
end