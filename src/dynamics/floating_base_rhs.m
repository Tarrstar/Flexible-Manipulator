function dx = floating_base_rhs( ...
        t, x, model, NB, dof, jointTorque)
%FLOATING_BASE_RHS State derivative for the floating-base manipulator.
% x = [quaternion(4); position(3); q; eta; spatialVelocity(6); qDot; etaDot].

positionSize = 7 + NB + dof;
theta = x(1:4);
theta = theta / norm(theta);
q = x(7 + (1:NB));
eta = x(7 + NB + (1:dof));
baseVelocity = x(positionSize + (1:6));
qDot = x(positionSize + 6 + (1:NB));
etaDot = x(positionSize + 6 + NB + (1:dof));

A = get_rotation_matrix_quaternion(theta);
baseState.velocity = baseVelocity;
baseState.rotation = A;
if isa(jointTorque, 'function_handle')
    appliedTorque = jointTorque(t);
else
    appliedTorque = jointTorque;
end
tau = [0; 0; 0; appliedTorque; appliedTorque; appliedTorque];
[qDDot, etaDDot, baseAcceleration] = FD_floating( ...
    model, q, qDot, tau, eta, etaDot, [], baseState);

[~, Gbar] = get_quaternion_G_matrices(theta);
thetaDot = 0.25 * Gbar.' * baseVelocity(1:3);
positionDot = A * baseVelocity(4:6);

dx = [thetaDot; positionDot; qDot; etaDot; ...
      baseAcceleration; qDDot; etaDDot{2}(:)];
end
