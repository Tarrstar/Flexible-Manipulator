function A = get_rotation_matrix(v, theta)
    symbolicAxis = is_casadi_type(v);
    symbolicAngle = is_casadi_type(theta);
    if numel(v) ~= 3 || (~symbolicAxis && ...
            (~isnumeric(v) || ~isreal(v) || ~isvector(v)))
        error('v has to be a 3d vector(3×1 or 1×3)');
    end

    if numel(theta) ~= 1 || (~symbolicAngle && ...
            (~isnumeric(theta) || ~isreal(theta) || ~isfinite(theta)))
        error('theta has to be a finite real unmber, its unit is rad');
    end

    if ~symbolicAxis && any(~isfinite(v))
        error('v cannot contain Inf or NaN');
    end

    v = v(:);

    norm_v = norm(v);

    if ~symbolicAxis && norm_v < eps
        error('v cannot be a zero vector');
    end

    v = v / norm_v;

    tilde_v = get_skew_sym_matrix(v);

    A = eye(3) ...
        + sin(theta) * tilde_v ...
        + (1 - cos(theta)) * (tilde_v * tilde_v);

end
