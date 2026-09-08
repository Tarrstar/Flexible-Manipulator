function X = local_deformation_xform(Phi, eta)

if isempty(Phi)
    X = eye(6);
    return;
end

[Ractive, ~, delta] = get_deformation_transform(Phi, eta);

% Ractive describes the physical rotation of the deformed endpoint.
% Xup requires the passive parent-to-child coordinate transformation.
E = Ractive.';
u = delta(4:6);

X = [E, zeros(3);
    -E * get_skew_sym_matrix(u), E];
end