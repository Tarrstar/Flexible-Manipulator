function [G_dot, G_bar_dot] = get_quaternion_G_dot_matrices(theta_dot)

    if ~isnumeric(theta_dot) || ~isreal(theta_dot) || ...
            ~isvector(theta_dot) || numel(theta_dot) ~= 4
        error('theta_dot must be a real 4D vector.');
    end

    if any(~isfinite(theta_dot))
        error('theta_dot cannot contain Inf or NaN.');
    end

    theta_dot = theta_dot(:);

    theta0_dot = theta_dot(1);
    theta1_dot = theta_dot(2);
    theta2_dot = theta_dot(3);
    theta3_dot = theta_dot(4);

    E_dot = [ ...
        -theta1_dot,  theta0_dot, -theta3_dot,  theta2_dot;
        -theta2_dot,  theta3_dot,  theta0_dot, -theta1_dot;
        -theta3_dot, -theta2_dot,  theta1_dot,  theta0_dot];

    E_bar_dot = [ ...
        -theta1_dot,  theta0_dot,  theta3_dot, -theta2_dot;
        -theta2_dot, -theta3_dot,  theta0_dot,  theta1_dot;
        -theta3_dot,  theta2_dot, -theta1_dot,  theta0_dot];

    G_dot     = 2 * E_dot;
    G_bar_dot = 2 * E_bar_dot;

end