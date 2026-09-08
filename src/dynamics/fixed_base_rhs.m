function dx = fixed_base_rhs( ...
        t, x, model, NB, dof, flexLink, forceInput)
%FIXED_BASE_RHS State derivative for a fixed-base flexible multibody model.
% x = [q; eta; qDot; etaDot].

positionSize = NB + dof;
q = x(1:NB);
eta = x(NB + (1:dof));
qDot = x(positionSize + (1:NB));
etaDot = x(positionSize + NB + (1:dof));

if isa(forceInput, 'function_handle')
    tau = forceInput(t);
else
    tau = forceInput;
end
tau = tau(:);
if numel(tau) ~= NB
    error('forceInput must return one torque per model link.');
end

[qDDot, etaDDot] = FD_flex( ...
    model, q, qDot, tau, eta, etaDot);
dx = [qDot; etaDot; qDDot; etaDDot{flexLink}(:)];
end
