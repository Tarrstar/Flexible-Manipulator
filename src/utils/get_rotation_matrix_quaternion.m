function A = get_rotation_matrix_quaternion(q)

    symbolic = is_casadi_type(q);
    if numel(q) ~= 4 || (~symbolic && (~isnumeric(q) || ~isreal(q) || ~isvector(q)))
        error('q has to be a real 4D quaternion vector (4x1 or 1x4).');
    end

    if ~symbolic && any(~isfinite(q))
        error('q cannot contain Inf or NaN.');
    end

    q = q(:);

    norm_q = norm(q);

    if ~symbolic && norm_q < eps
        error('q cannot be a zero quaternion.');
    end

    % Normalize, so non-unit inputs are also accepted.
    q = q / norm_q;

    qw = q(1);
    qx = q(2);
    qy = q(3);
    qz = q(4);

    A = [ ...
        1 - 2*(qy^2 + qz^2),     2*(qx*qy - qw*qz),     2*(qx*qz + qw*qy);
            2*(qx*qy + qw*qz), 1 - 2*(qx^2 + qz^2),     2*(qy*qz - qw*qx);
            2*(qx*qz - qw*qy),     2*(qy*qz + qw*qx), 1 - 2*(qx^2 + qy^2)];
end
