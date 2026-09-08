function [Irr, Ire, Mee, C, c] = floating_frame_terms(body, v, eta, eta_d)
% Convert FlexBody's [quaternion; translation; eta] equations to spatial
% floating-frame blocks. v is ordered as [omega; linearVelocity] in the
% floating frame.

    qBody = body.q(:);
    v = v(:);
    eta = eta(:);
    eta_d = eta_d(:);

    qBody(8:end) = eta;

    theta = qBody(1:4);
    [~, Gbar] = get_quaternion_G_matrices(theta);
    A = get_rotation_matrix_quaternion(theta);

    omega = v(1:3);
    linearVelocity = v(4:6);

    % Gbar = 2*Ebar, therefore theta_dot = 0.25*Gbar'*omega.
    thetaDot = 0.25 * Gbar.' * omega;
    qdotBody = [thetaDot;
                A * linearVelocity;
                eta_d];

    M = body.get_mass_matrix(qBody);
    Qv = body.get_quadratic_velocity_vec(qBody, qdotBody);
    Qe = body.get_elastic_force(qBody, qdotBody);

    Jr = [0.25 * Gbar.', zeros(4, 3);
          zeros(3, 3), A];

    % Spatial acceleration already includes the translational convective
    % term. The remaining coordinate-map bias is quaternion-only.
    [~, GbarDot] = get_quaternion_G_dot_matrices(thetaDot);
    qddBias = [
    0.25 * GbarDot.' * omega;
    A * cross(omega, linearVelocity);
    zeros(body.nEta, 1)
];
    bias = -Qv + M * qddBias;

    Mrr = M(1:7, 1:7);
    Irr = Jr.' * Mrr * Jr;

    nEta = body.nEta;
    if nEta == 0
        Ire = zeros(6, 0);
        Mee = zeros(0, 0);
    else
        Ire = reshape(Jr.' * M(1:7, 8:(7 + nEta)), 6, nEta);
        Mee = reshape(M(8:(7 + nEta), 8:(7 + nEta)), nEta, nEta);
    end

    C = Jr.' * bias(1:7);
    c = bias(8:end) + Qe(8:end);
end
