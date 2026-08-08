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

%% ================== 1. Beam parameters ==================

Lbeam  = 2.0;
W_beam = 0.1;
T_beam = 0.05;

rho = 2700.0;
E   = 7.0e7;
nu  = 0.33;

nElements = 2;


beam.length = Lbeam;
beam.width = W_beam;
beam.thickness = T_beam;
beam.density = rho;
beam.youngs_modulus = E;
beam.poisson_ratio = nu;
beam.elements_num = nElements;
beam.fixed_root = true;

femBeam = FEM_Beam(beam);

[~, K_beam] = femBeam.get_matrices_3d();

[feShapeMatrix, stiffnessDofs] = femBeam.build_direct_fe_shape_matrix();

nFlexibleDofs = size(feShapeMatrix, 2);
Kff = full(K_beam(stiffnessDofs, stiffnessDofs));
Kff = 0.5 * (Kff + Kff.');

if size(Kff, 1) ~= nFlexibleDofs
    error('Kff size does not match the direct FE coordinate size.');
end


% Numerical quadrature for the FFRF mass integrals.
%nx = max(2 * nElements, 24);
nx = 4;
ny = 3;
nz = 3;

[u0, S, weights] = femBeam.build_ffrf_data( ...
    feShapeMatrix, nx, ny, nz);

% No structural damping, same as Exudyn example.
Dff = zeros(nFlexibleDofs);

params.density = rho;
params.u0      = u0;
params.S       = S;
params.weights = weights;
params.Kff     = Kff;
params.Dff     = Dff;

% q = [theta; R; eta]
params.q = [1; 0; 0; 0; ...
            0; 0; 0; ...
            zeros(nFlexibleDofs, 1)];

params.qdot = zeros(7 + nFlexibleDofs, 1);

body = FlexBody(params);

nq = 7 + nFlexibleDofs;

%% ================== 4. Simulation settings ==================

simulationTime = 2.0;
timeStep = 2e-3;

tVec = 0:timeStep:simulationTime;

% State:
%
% x = [q;
%      qdot]
%
x0 = [params.q;
      params.qdot];

% Applied torque at the left reference point.
% Same as:
%
% LoadTorqueVector(..., loadVector=[0, 0, 1.0])
%
Mz = 1.0;

%% ================== 5. Time integration ==================

odeOptions = odeset( ...
    'RelTol', 1e-6, ...
    'AbsTol', 1e-8, ...
    'MaxStep', timeStep);

rhs = @(t, x) rhs_flexbody(t, x, body, Mz);
[~, xHist] = ode45(rhs, tVec, x0, odeOptions);

log_data = zeros(numel(tVec), 6);

for k = 1:numel(tVec)
    q = xHist(k, 1:nq).';
    qdot = xHist(k, nq+1:end).';

    [q, qdot] = normalize_quaternion_state(q, qdot);
    xHist(k, 1:nq) = q.';
    xHist(k, nq+1:end) = qdot.';

    [base_x, base_y, base_theta, tip_x, tip_y] = ...
        get_output_for_plot(q, femBeam, feShapeMatrix);

    log_data(k, :) = [tVec(k), base_x, base_y, ...
                      base_theta, tip_x, tip_y];
end

%% ================== 6. Save and plot ==================

csvName = fullfile(dataDir, 'floating_beam_flexbody_fe_ode45.csv');

writematrix_with_header(csvName, log_data, ...
    'time,base_x,base_y,base_theta,tip_x,tip_y');

exudynCsvName = fullfile(dataDir, 'floating_beam_exudyn.csv');
rigidCsvName = fullfile(dataDir, 'floating_rigid_beam_exudyn.csv');
comparisonFigName = fullfile(dataDir, ...
    'floating_beam_exudyn_vs_flexbody.png');

exudyn_data = [];
rigid_data = [];

if exist(exudynCsvName, 'file')
    exudyn_data = readmatrix(exudynCsvName);
    exudyn_data = exudyn_data(all(isfinite(exudyn_data), 2), :);

else
    warning(['Exudyn CSV was not found. Run refs/floating_beam.py ', ...
             'first to create: %s'], exudynCsvName);
end

if exist(rigidCsvName, 'file')
    rigid_data = readmatrix(rigidCsvName);
    rigid_data = rigid_data(all(isfinite(rigid_data), 2), :);

else
    warning(['Rigid-body CSV was not found. Run refs/floating_rigid_body.py ', ...
             'first to create: %s'], rigidCsvName);
end

figure('Color', 'w', 'Position', [120, 100, 1050, 720]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(log_data(:, 1), log_data(:, 2), 'LineWidth', 1.2); hold on;
plot(log_data(:, 1), log_data(:, 3), 'LineWidth', 1.2);
baseLegend = {'MATLAB Base X', 'MATLAB Base Y'};


grid on;
title('Base');
xlabel('Time / s');
ylabel('Base position');
legend(baseLegend{:}, 'Location', 'best');

nexttile;
plot(log_data(:, 1), log_data(:, 5), 'LineWidth', 1.2); hold on;
plot(log_data(:, 1), log_data(:, 6), 'LineWidth', 1.2);
tipLegend = {'MATLAB Tip X', 'MATLAB Tip Y'};


grid on;
title('Tip');
xlabel('Time / s');
ylabel('Tip position');
legend(tipLegend{:}, 'Location', 'best');

saveas(gcf, comparisonFigName);
disp(['Comparison figure saved to ', comparisonFigName]);

disp(['Simulation finished. Data saved to ', csvName]);

%% ========================================================================
% Local functions
% ========================================================================

function [base_x, base_y, base_theta, tip_x, tip_y] = ...
        get_output_for_plot(q, femBeam, feShapeMatrix)

    theta = q(1:4);
    R = q(5:7);
    eta = q(8:end);

    A = get_rotation_matrix_quaternion(theta);

    base_x = R(1);
    base_y = R(2);
    base_theta = atan2(A(2, 1), A(1, 1));

    tipFeDisp = femBeam.evaluate_displacement( ...
        feShapeMatrix, femBeam.length, 0, 0);

    rTipLocal = [femBeam.length; 0; 0] + tipFeDisp * eta;
    rTipGlobal = R + A * rTipLocal;

    tip_x = rTipGlobal(1);
    tip_y = rTipGlobal(2);
end
