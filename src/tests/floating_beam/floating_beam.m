clear; clc; close all;

%% ================== 1. Beam parameters ==================

Lbeam  = 2.0;
W_beam = 0.1;
T_beam = 0.05;

rho = 2700.0;
E   = 7.0e7;

A_cs = W_beam * T_beam;
Izz  = W_beam * T_beam^3 / 12.0;

EA = E * A_cs;
EI = E * Izz;

% Same idea as Exudyn example:
% Cable2D uses a planar beam. Here we keep only planar bending
% in the local x-y plane.
nModes = 5;

betaL_all = [ ...
    1.875104068711961, ...
    4.694091132974174, ...
    7.854757438237612, ...
    10.99554073487547, ...
    14.13716839104647, ...
    17.27875953208824, ...
    20.42035225104125];

betaL = betaL_all(1:nModes);

%% ================== 2. Build AMM shape data ==================

% Numerical quadrature for the FFRF mass integrals.
nx = 16;
ny = 3;
nz = 3;

phiFcn = @(x) local_phi(x, Lbeam, betaL);
dPhiFcn = @(x) local_dphi(x, Lbeam, betaL);

[u0, S, weights] = build_ffrf_beam_bending_y_data( ...
    Lbeam, W_beam, T_beam, nx, ny, nz, phiFcn, dPhiFcn);

% Modal stiffness from Euler-Bernoulli bending energy:
%
% Kff = integral_0^L EI * Phi''(x)^T * Phi''(x) dx
%
Kff = build_bending_stiffness_matrix(Lbeam, EI, betaL, 40);

% No structural damping, same as Exudyn example.
Dff = zeros(nModes);

%% ================== 3. Create FlexBody ==================

params.density = rho;
params.u0      = u0;
params.S       = S;
params.weights = weights;
params.Kff     = Kff;
params.Dff     = Dff;

% q = [theta; R; eta]
params.q = [1; 0; 0; 0; ...
            0; 0; 0; ...
            zeros(nModes, 1)];

params.qdot = zeros(7 + nModes, 1);

body = FlexBody(params);

nq = 7 + nModes;

%% ================== 4. Simulation settings ==================

simulationTime = 2.0;
timeStep = 2e-3;

tVec = 0:timeStep:simulationTime;
nSteps = numel(tVec);

% State:
%
% x = [q;
%      qdot]
%
x = [params.q;
     params.qdot];

log_data = zeros(nSteps, 6);

% Applied torque at the left reference point.
% Same as:
%
% LoadTorqueVector(..., loadVector=[0, 0, 1.0])
%
Mz = 1.0;

%% ================== 5. Time integration ==================

for k = 1:nSteps
    t = tVec(k);

    q    = x(1:nq);
    qdot = x(nq+1:end);

    [q, qdot] = normalize_quaternion_state(q, qdot);
    x(1:nq) = q;
    x(nq+1:end) = qdot;

    [base_x, base_y, base_theta, tip_x, tip_y] = ...
        get_output_for_plot(body, q, Lbeam, betaL);

    log_data(k, :) = [t, base_x, base_y, base_theta, tip_x, tip_y];

    if k < nSteps
        % Fixed-step RK4, close to Exudyn fixed time step output.
        h = timeStep;

        k1 = rhs_flexbody(t, x, body, Mz);
        k2 = rhs_flexbody(t + 0.5*h, x + 0.5*h*k1, body, Mz);
        k3 = rhs_flexbody(t + 0.5*h, x + 0.5*h*k2, body, Mz);
        k4 = rhs_flexbody(t + h, x + h*k3, body, Mz);

        x = x + h/6 * (k1 + 2*k2 + 2*k3 + k4);

        q    = x(1:nq);
        qdot = x(nq+1:end);

        [q, qdot] = normalize_quaternion_state(q, qdot);

        x(1:nq) = q;
        x(nq+1:end) = qdot;
    end
end

%% ================== 6. Save and plot ==================

csvName = 'floating_beam_flexbody_amm.csv';

writematrix_with_header(csvName, log_data, ...
    'time,base_x,base_y,base_theta,tip_x,tip_y');

figure('Color', 'w', 'Position', [100, 100, 900, 700]);

subplot(3, 1, 1);
plot(log_data(:, 1), log_data(:, 2), 'LineWidth', 1.2); hold on;
plot(log_data(:, 1), log_data(:, 3), 'LineWidth', 1.2);
grid on;
legend('Base X', 'Base Y');
ylabel('Base position');

subplot(3, 1, 2);
plot(log_data(:, 1), log_data(:, 4), 'g', 'LineWidth', 1.2);
grid on;
legend('Base Theta');
ylabel('\theta / rad');

subplot(3, 1, 3);
plot(log_data(:, 1), log_data(:, 5), 'LineWidth', 1.2); hold on;
plot(log_data(:, 1), log_data(:, 6), 'LineWidth', 1.2);
grid on;
legend('Tip X', 'Tip Y');
xlabel('Time / s');
ylabel('Tip position');

disp(['Simulation finished. Data saved to ', csvName]);

%% ========================================================================
% Local functions
% ========================================================================

function dx = rhs_flexbody(~, x, body, Mz)

    nq = numel(x) / 2;

    q    = x(1:nq);
    qdot = x(nq+1:end);

    [q, qdot] = normalize_quaternion_state(q, qdot);

    M  = body.get_mass_matrix(q);
    Qv = body.get_quadratic_velocity_vec(q, qdot);
    Qe = body.get_elastic_force(q, qdot);

    theta = q(1:4);
    thetaDot = qdot(1:4);

    [G, G_bar] = get_quaternion_G_matrices(theta);

    % External generalized force.
    %
    % The physical torque is about the global z axis.
    % For a global torque n, power is:
    %
    %   P = n^T * omega = n^T * G * thetaDot
    %
    % Therefore:
    %
    %   Q_theta = G^T * n
    %
    Qext = zeros(nq, 1);

    torque_global = [0; 0; Mz];
    Qext(1:4) = G.' * torque_global;

    % No external force on the reference point.
    Qext(5:7) = [0; 0; 0];

    % Quaternion acceleration constraint:
    %
    % theta^T theta = 1
    % 2 theta^T thetaDDot + 2 thetaDot^T thetaDot = 0
    %
    Phi_q = [2 * theta.', zeros(1, nq - 4)];
    gamma = -2 * (thetaDot.' * thetaDot);

    % Dynamics:
    %
    % M qDDot + Qv + Qe = Qext + Phi_q^T lambda
    %
    % Therefore:
    %
    % [M, -Phi_q^T; Phi_q, 0] [qDDot; lambda]
    % =
    % [Qext - Qv - Qe; gamma]
    %
    KKT = [M, -Phi_q.';
           Phi_q, 0];

    rhs = [Qext - Qv - Qe;
           gamma];

    sol = KKT \ rhs;

    qddot = sol(1:nq);

    dx = [qdot;
          qddot];
end


function [base_x, base_y, base_theta, tip_x, tip_y] = ...
        get_output_for_plot(~, q, Lbeam, betaL)

    theta = q(1:4);
    R = q(5:7);
    eta = q(8:end);

    A = get_rotation_matrix_quaternion(theta);

    % Base position is the reference origin.
    base_x = R(1);
    base_y = R(2);

    % For planar rotation about z:
    %
    % A = [cos(theta), -sin(theta), 0;
    %      sin(theta),  cos(theta), 0;
    %      0,           0,          1]
    %
    base_theta = atan2(A(2, 1), A(1, 1));

    % Tip centerline local position:
    %
    % x = L
    % y = Phi(L)*eta
    % z = 0
    %
    PhiL = local_phi(Lbeam, Lbeam, betaL);

    r_tip_local = [Lbeam;
                   PhiL * eta;
                   0];

    r_tip_global = R + A * r_tip_local;

    tip_x = r_tip_global(1);
    tip_y = r_tip_global(2);
end

function [u0, S, weights] = build_ffrf_beam_bending_y_data( ...
        L, W, T, nx, ny, nz, phiFcn, dPhiFcn)

    [xi_x, wi_x] = gauss_legendre_1d(nx);
    [xi_y, wi_y] = gauss_legendre_1d(ny);
    [xi_z, wi_z] = gauss_legendre_1d(nz);

    xNodes = L / 2 * (xi_x + 1);
    yNodes = T / 2 * xi_y;
    zNodes = W / 2 * xi_z;

    jac = L / 2 * T / 2 * W / 2;

    PhiTest = phiFcn(xNodes(1));
    PhiTest = PhiTest(:).';

    nEta = numel(PhiTest);
    nGauss = nx * ny * nz;

    u0 = zeros(3, nGauss);
    S = zeros(3, nEta, nGauss);
    weights = zeros(nGauss, 1);

    g = 0;

    for ix = 1:nx
        x = xNodes(ix);

        Phi = phiFcn(x);
        Phi = Phi(:).';

        dPhi = dPhiFcn(x);
        dPhi = dPhi(:).';

        for iy = 1:ny
            y = yNodes(iy);

            for iz = 1:nz
                z = zNodes(iz);

                g = g + 1;

                u0(:, g) = [x; y; z];

                S(:, :, g) = [ ...
                    -y * dPhi;
                     Phi;
                     zeros(1, nEta)];

                weights(g) = wi_x(ix) * wi_y(iy) * wi_z(iz) * jac;
            end
        end
    end
end

