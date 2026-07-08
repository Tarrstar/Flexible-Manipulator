function [M_beam, K_beam] = get_beam_matrices_3d(beam)


    L = beam.length;
    n_element = beam.elements_num;
    l = L / n_element;

    width = beam.width;
    thickness = beam.thickness;

    rho = beam.density;
    E = beam.youngs_modulus;

    if isfield(beam, 'poisson_ratio')
        nu = beam.poisson_ratio;
    else
        nu = 0.33;
    end

    if isfield(beam, 'shear_modulus')
        G = beam.shear_modulus;
    else
        G = E / (2 * (1 + nu));
    end

    if isfield(beam, 'fixed_root')
        fixed_root = beam.fixed_root;
    else
        fixed_root = false;
    end

    % Section properties

    A = width * thickness;

    % Local coordinate:
    % x : beam axis
    % y : width direction
    % z : thickness direction
    %
    % Iyy is used for bending in x-z plane, i.e. w-theta_y.
    % Izz is used for bending in x-y plane, i.e. v-theta_z.
    Iyy = width * thickness^3 / 12;
    Izz = thickness * width^3 / 12;

    % Approximate Saint-Venant torsion constant for rectangular section.
    J = get_rectangular_torsion_constant(width, thickness);

    % Polar second moment used for torsional consistent mass.
    Ip = Iyy + Izz;

    % Element stiffness and mass matrices

    Ke = zeros(12, 12);
    Me = zeros(12, 12);

    % Element DOF order:
    %
    % [u1, v1, w1, theta_x1, theta_y1, theta_z1, ...
    %  u2, v2, w2, theta_x2, theta_y2, theta_z2]

    % 1. Axial part: u

    Ke_u = E * A / l * ...
        [ 1, -1;
         -1,  1];

    Me_u = rho * A * l / 6 * ...
        [2, 1;
         1, 2];

    idx_u = [1, 7];

    Ke(idx_u, idx_u) = Ke(idx_u, idx_u) + Ke_u;
    Me(idx_u, idx_u) = Me(idx_u, idx_u) + Me_u;

    % 2. Torsional part: theta_x

    Ke_tx = G * J / l * ...
        [ 1, -1;
         -1,  1];

    Me_tx = rho * Ip * l / 6 * ...
        [2, 1;
         1, 2];

    idx_tx = [4, 10];

    Ke(idx_tx, idx_tx) = Ke(idx_tx, idx_tx) + Ke_tx;
    Me(idx_tx, idx_tx) = Me(idx_tx, idx_tx) + Me_tx;

    % 3. Bending in x-y plane: v and theta_z

    Ke_v = E * Izz / l^3 * ...
        [ 12,    6*l,   -12,    6*l;
          6*l,  4*l^2, -6*l,  2*l^2;
         -12,   -6*l,    12,   -6*l;
          6*l,  2*l^2, -6*l,  4*l^2];

    Me_v = rho * A * l / 420 * ...
        [ 156,    22*l,    54,   -13*l;
          22*l,  4*l^2,  13*l,  -3*l^2;
           54,    13*l,   156,   -22*l;
         -13*l, -3*l^2, -22*l,   4*l^2];

    idx_v = [2, 6, 8, 12];

    Ke(idx_v, idx_v) = Ke(idx_v, idx_v) + Ke_v;
    Me(idx_v, idx_v) = Me(idx_v, idx_v) + Me_v;

    % 4. Bending in x-z plane: w and theta_y

    Ke_w_standard = E * Iyy / l^3 * ...
        [ 12,    6*l,   -12,    6*l;
          6*l,  4*l^2, -6*l,  2*l^2;
         -12,   -6*l,    12,   -6*l;
          6*l,  2*l^2, -6*l,  4*l^2];

    Me_w_standard = rho * A * l / 420 * ...
        [ 156,    22*l,    54,   -13*l;
          22*l,  4*l^2,  13*l,  -3*l^2;
           54,    13*l,   156,   -22*l;
         -13*l, -3*l^2, -22*l,   4*l^2];

    % Standard Hermite matrix uses slope beta = dw/dx.
    % But physical small rotation is theta_y = -dw/dx.
    % Therefore:
    %
    %   beta = -theta_y
    %
    Tsign = diag([1, -1, 1, -1]);

    Ke_w = Tsign' * Ke_w_standard * Tsign;
    Me_w = Tsign' * Me_w_standard * Tsign;

    idx_w = [3, 5, 9, 11];

    Ke(idx_w, idx_w) = Ke(idx_w, idx_w) + Ke_w;
    Me(idx_w, idx_w) = Me(idx_w, idx_w) + Me_w;

    % Assemble global sparse matrices

    n_node = n_element + 1;
    dof_per_node = 6;
    total_dof = dof_per_node * n_node;

    I_all = [];
    J_all = [];
    K_all = [];
    M_all = [];

    for e = 1:n_element

        node1 = e;
        node2 = e + 1;

        dofs = [ ...
            get_node_dofs(node1), ...
            get_node_dofs(node2)];

        [ii, jj] = ndgrid(dofs, dofs);

        I_all = [I_all; ii(:)];
        J_all = [J_all; jj(:)];

        K_all = [K_all; Ke(:)];
        M_all = [M_all; Me(:)];
    end

    K_beam = sparse(I_all, J_all, K_all, total_dof, total_dof);
    M_beam = sparse(I_all, J_all, M_all, total_dof, total_dof);

    K_beam = 0.5 * (K_beam + K_beam');
    M_beam = 0.5 * (M_beam + M_beam');

    % Remove fixed root DOFs if needed

    if fixed_root
        fixed_dofs = 1:6;
        free_dofs = setdiff(1:total_dof, fixed_dofs);

        K_beam = K_beam(free_dofs, free_dofs);
        M_beam = M_beam(free_dofs, free_dofs);
    end

end

% ========================================================================
% Local functions
% ========================================================================

function dofs = get_node_dofs(node_id)

    dof_per_node = 6;

    start_id = dof_per_node * (node_id - 1) + 1;

    dofs = start_id:(start_id + dof_per_node - 1);

end

function J = get_rectangular_torsion_constant(width, thickness)
% Approximate Saint-Venant torsion constant of a rectangular section.

    a = max(width, thickness);
    b = min(width, thickness);

    J = a * b^3 * ...
        (1/3 - 0.21 * (b/a) * (1 - b^4 / (12*a^4)));

end