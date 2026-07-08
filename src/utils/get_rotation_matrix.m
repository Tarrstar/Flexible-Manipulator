function A = get_rotation_matrix(v, theta)
    if ~isnumeric(v) || ~isreal(v) || ~isvector(v) || numel(v) ~= 3
        error('v has to be a 3d vector(3×1 or 1×3)');
    end

    if ~isnumeric(theta) || ~isreal(theta) || ~isscalar(theta) || ...
            ~isfinite(theta)
        error('theta has to be a finite real unmber, its unit is rad');
    end

    if any(~isfinite(v))
        error('v cannot contain Inf or NaN');
    end

    v = v(:);

    norm_v = norm(v);

    if norm_v < eps
        error('v cannot be a zero vector');
    end

    v = v / norm_v;

    tilde_v = get_skew_sym_matrix(v);

    A = eye(3) ...
        + sin(theta) * tilde_v ...
        + (1 - cos(theta)) * (tilde_v * tilde_v);

end