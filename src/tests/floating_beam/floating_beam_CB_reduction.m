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
nInteriorModes = 6;

beam.length = Lbeam;
beam.width = W_beam;
beam.thickness = T_beam;
beam.density = rho;
beam.youngs_modulus = E;
beam.poisson_ratio = nu;
beam.elements_num = nElements;
beam.fixed_root = true;
beam.nInteriorModes = nInteriorModes;
beam.boundaryNodeIds = [1, nElements + 1];

femBeam = FEM_Beam_CB_Reduction(beam);

%% ================== 2. Craig-Bampton coordinates ==================

nx = 4;
ny = 3;
nz = 3;

[u0, S, weights, cbShapeMatrix, Kff, reductionMatrix, cbInfo] = ...
    femBeam.build_ffrf_data(nInteriorModes, nx, ny, nz);

nFlexibleDofs = size(cbShapeMatrix, 2);
Dff = zeros(nFlexibleDofs);

disp(['CB boundary DOFs: ', num2str(cbInfo.nBoundaryDofs)]);
disp(['CB interior modes: ', num2str(cbInfo.nInteriorModes)]);
disp(['CB reduced coordinates: ', num2str(size(reductionMatrix, 2))]);

%% ================== 3. Create FlexBody ==================

params.density = rho;
params.u0      = u0;
params.S       = S;
params.weights = weights;
params.Kff     = Kff;
params.Dff     = Dff;

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
x0 = [params.q;
      params.qdot];

F = [0;0;0;0;1;2];

%% ================== 5. Time integration ==================

odeOptions = odeset( ...
    'RelTol', 1e-6, ...
    'AbsTol', 1e-8, ...
    'MaxStep', timeStep);

rhs = @(t, x) rhs_flexbody(t, x, body, F);
[~, xHist] = ode45(rhs, tVec, x0, odeOptions);

log_data = zeros(numel(tVec), 8);

for k = 1:numel(tVec)
    q = xHist(k, 1:nq).';
    qdot = xHist(k, nq+1:end).';

    [q, qdot] = normalize_quaternion_state(q, qdot);
    xHist(k, 1:nq) = q.';
    xHist(k, nq+1:end) = qdot.';

    [base_x, base_y, base_z, base_theta, tip_x, tip_y, tip_z] = ...
        get_output_for_plot(q, femBeam, cbShapeMatrix);

    log_data(k, :) = [tVec(k), base_x, base_y, base_z, ...
                      base_theta, tip_x, tip_y, tip_z];
end

%% ================== 6. Save and compare ==================

csvName = fullfile(dataDir, 'floating_beam_flexbody_cb_ode45.csv');

writematrix_with_header(csvName, log_data, ...
    'time,base_x,base_y,base_z,base_theta,tip_x,tip_y,tip_z');

exudynCsvName = fullfile(dataDir, 'floating_beam_3d_exudyn.csv');
comparisonFigName = fullfile(dataDir, ...
    'floating_beam_cb_3d_exudyn_vs_flexbody.png');

exudyn_data = [];

if exist(exudynCsvName, 'file')
    exudyn_data = readmatrix(exudynCsvName);
    exudyn_data = exudyn_data(all(isfinite(exudyn_data), 2), :);

    if size(exudyn_data, 2) < 8
        warning('Exudyn CSV does not contain the expected 8 columns: %s', ...
            exudynCsvName);
        exudyn_data = [];
    end
else
    warning(['Exudyn CSV was not found. Run refs/floating_beam_3d.py ', ...
             'first to create: %s'], exudynCsvName);
end

figure('Color', 'w', 'Position', [120, 100, 1050, 720]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(log_data(:, 1), log_data(:, 2), 'LineWidth', 1.2); hold on;
plot(log_data(:, 1), log_data(:, 3), 'LineWidth', 1.2);
plot(log_data(:, 1), log_data(:, 4), 'LineWidth', 1.2);
baseLegend = {'CB Base X', 'CB Base Y', 'CB Base Z'};

if ~isempty(exudyn_data)
    plot(exudyn_data(:, 1), exudyn_data(:, 2), '--', 'LineWidth', 1.2);
    plot(exudyn_data(:, 1), exudyn_data(:, 3), '--', 'LineWidth', 1.2);
    plot(exudyn_data(:, 1), exudyn_data(:, 4), '--', 'LineWidth', 1.2);
    baseLegend = [baseLegend, ...
        {'Exudyn Base X', 'Exudyn Base Y', 'Exudyn Base Z'}];
end

grid on;
% title('Base');
xlabel('Time (s)');
ylabel('Base position');
legend(baseLegend{:}, 'Location', 'best');

nexttile;
plot(log_data(:, 1), log_data(:, 6), 'LineWidth', 1.2); hold on;
plot(log_data(:, 1), log_data(:, 7), 'LineWidth', 1.2);
plot(log_data(:, 1), log_data(:, 8), 'LineWidth', 1.2);
tipLegend = {'CB Tip X', 'CB Tip Y', 'CB Tip Z'};

if ~isempty(exudyn_data)
    plot(exudyn_data(:, 1), exudyn_data(:, 6), '--', 'LineWidth', 1.2);
    plot(exudyn_data(:, 1), exudyn_data(:, 7), '--', 'LineWidth', 1.2);
    plot(exudyn_data(:, 1), exudyn_data(:, 8), '--', 'LineWidth', 1.2);
    tipLegend = [tipLegend, ...
        {'Exudyn Tip X', 'Exudyn Tip Y', 'Exudyn Tip Z'}];
end

grid on;
% title('Tip');
xlabel('Time (s)');
ylabel('Tip position');
legend(tipLegend{:}, 'Location', 'best');

saveas(gcf, comparisonFigName);

disp(['CB FlexBody data saved to ', csvName]);
disp(['CB comparison figure saved to ', comparisonFigName]);

%% ========================================================================
% Local functions
% ========================================================================

function [base_x, base_y, base_z, base_theta, tip_x, tip_y, tip_z] = ...
        get_output_for_plot(q, femBeam, cbShapeMatrix)

    theta = q(1:4);
    R = q(5:7);
    eta = q(8:end);

    A = get_rotation_matrix_quaternion(theta);

    base_x = R(1);
    base_y = R(2);
    base_z = R(3);
    base_theta = atan2(A(2, 1), A(1, 1));

    tipFeDisp = femBeam.evaluate_displacement( ...
        cbShapeMatrix, femBeam.length, 0, 0);

    rTipLocal = [femBeam.length; 0; 0] + tipFeDisp * eta;
    rTipGlobal = R + A * rTipLocal;

    tip_x = rTipGlobal(1);
    tip_y = rTipGlobal(2);
    tip_z = rTipGlobal(3);
end
