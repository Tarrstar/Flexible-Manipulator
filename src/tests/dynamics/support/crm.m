function X = crm(v)
% Spatial motion cross-product operator used by the dynamics test.
v = v(:);
omegaTilde = get_skew_sym_matrix(v(1:3));
linearTilde = get_skew_sym_matrix(v(4:6));
X = [omegaTilde, zeros(3); linearTilde, omegaTilde];
end
