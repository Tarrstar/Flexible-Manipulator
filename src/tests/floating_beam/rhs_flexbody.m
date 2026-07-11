function dx = rhs_flexbody(~, x, body, F)

    nq = numel(x) / 2;

    q = x(1:nq);
    qdot = x(nq+1:end);

    [q, qdot] = normalize_quaternion_state(q, qdot);

    M = body.get_mass_matrix(q);
    Qv = body.get_quadratic_velocity_vec(q, qdot);
    Qe = body.get_elastic_force(q, qdot);

    theta = q(1:4);
    thetaDot = qdot(1:4);

    [G, ~] = get_quaternion_G_matrices(theta);

    Qext = zeros(nq, 1);
    torque_global = [F(1); F(2); F(3)];
    Qext(1:4) = G.' * torque_global;
    Qext(5:7) = [F(4); F(5); F(6)];

    Phi_q = [2 * theta.', zeros(1, nq - 4)];
    gamma = -2 * (thetaDot.' * thetaDot);

    KKT = [M, -Phi_q.';
           Phi_q, 0];

    rhs = [Qext + Qv - Qe;
           gamma];

    sol = KKT \ rhs;
    qddot = sol(1:nq);

    dx = [qdot;
          qddot];
end
