function X = local_deformation_xform(Phi, eta)
if isempty(Phi)
    X = eye(6);
    return;
end
[R, ~, delta] = get_deformation_transform(Phi, eta);
u = delta(4:6);
X = [R, zeros(3); -R * get_skew_sym_matrix(u), R];
end