function X = xlt(r)
% Spatial translation transform used by jcalc for the test's Pz joint.
r = r(:);
S = get_skew_sym_matrix(r);
X = [eye(3), zeros(3); -S, eye(3)];
end
