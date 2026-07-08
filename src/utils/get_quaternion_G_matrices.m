function [G, G_bar] = get_quaternion_G_matrices(q)

    if ~isnumeric(q) || ~isreal(q) || ...
            ~isvector(q) || numel(q) ~= 4
        error(['q must be a real 4D Euler-parameter vector ', ...
               '(4x1 or 1x4).']);
    end

    if any(~isfinite(q))
        error('q cannot contain Inf or NaN.');
    end

    q = q(:);

    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);

    E = [ ...
        -q1,  q0, -q3,  q2;
        -q2,  q3,  q0, -q1;
        -q3, -q2,  q1,  q0];

    E_bar = [ ...
        -q1,  q0,  q3, -q2;
        -q2, -q3,  q0,  q1;
        -q3,  q2, -q1,  q0];

    G     = 2 * E;
    G_bar = 2 * E_bar;

end