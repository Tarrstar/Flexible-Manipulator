function [R, T, delta] = get_deformation_transform(Phi, eta)

    eta = eta(:);
    delta = Phi * eta;

    theta = delta(1:3);
    u     = delta(4:6);
    theta_tilde = get_skew_sym_matrix(theta);
    if is_casadi_type(theta)

        angleSquared = theta.' * theta;
        small = angleSquared < 1e-16;
        safeSquared = if_else(small, 1, angleSquared);
        safeAngle = sqrt(safeSquared);
        A = if_else(small, 1 - angleSquared/6 + angleSquared^2/120, ...
            sin(safeAngle)/safeAngle);
        B = if_else(small, 0.5 - angleSquared/24 + angleSquared^2/720, ...
            (1 - cos(safeAngle))/safeSquared);
    else
        angleSquared = theta.' * theta;
        if angleSquared < 1e-16
            A = 1-angleSquared/6+angleSquared^2/120;
            B = 0.5-angleSquared/24+angleSquared^2/720;
        else
            angle = sqrt(angleSquared);
            A = sin(angle) / angle;
            B = (1 - cos(angle)) / angleSquared;
        end
    end

    R = eye(3) ...
        + A * theta_tilde ...
        + B * theta_tilde^2;

    T = [R,       u;
         0, 0, 0, 1];

end
