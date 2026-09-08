function [qdd, eta_dd, base_accel] = FD_floating( ...
        model, q, qd, tau, eta, eta_d, f_ext, base_state)
%FD_FLOATING Forward dynamics for a single-root floating flexible tree.
% Link 1 is always the floating base. Its q, qd, and tau entries are
% placeholders, and the third forward recursion starts with link 2.

NB = model.NB;

baseTwist = local_base_field(base_state, 'velocity', zeros(6, 1));
baseWrench = local_base_field(base_state, 'wrench', zeros(6, 1));

q = q(:); qd = qd(:); tau = tau(:);

fixedJoint = false(NB, 1);
if isfield(model, 'fixed_joints') && ~isempty(model.fixed_joints)
    fixedSpec = model.fixed_joints;
    if islogical(fixedSpec) && numel(fixedSpec) == NB
        fixedJoint = fixedSpec(:);
    else
        fixedIds = fixedSpec(:);
        fixedJoint(fixedIds) = true;
    end
end
fixedJoint(1) = true;

if nargin < 7 || isempty(f_ext)
    f_ext = cell(NB, 1);
elseif ~iscell(f_ext)
    f_ext = {f_ext};
elseif numel(f_ext) < NB
    f_ext(end + 1:NB) = {[]};
end

if isfield(model, 'gravity')
    a_grav = model.gravity(:);
else
    a_grav = zeros(6, 1);
end
if numel(a_grav) == 3
    a_grav = [0; 0; 0; a_grav];
end

if isfield(base_state, 'rotation') && ~isempty(base_state.rotation)
    A0 = base_state.rotation;
    a_grav = [A0.' * a_grav(1:3); A0.' * a_grav(4:6)];
end

% Forward pass (Pass 1: Kinematics and local dynamics terms)
for i = 1:NB
    if is_flexible_link(model, i)
        [B{i}, eta_i, eta_d_i] = local_modal_state(model, i, eta, eta_d);
    else
        B{i} = zeros(6, 0);
        eta_i = zeros(0, 1); eta_d_i = zeros(0, 1);
    end
    etaState{i} = eta_i;
    etaDotState{i} = eta_d_i;

    if i == 1
        Xup{i} = eye(6);
        S{i} = zeros(6, 1);
        edgeR{i} = eye(6);
        edgeB{i} = zeros(6, 0);
        v{i} = baseTwist;
        gamma{i} = zeros(6, 1);
    else
        if fixedJoint(i)
            [XJ, ~] = jcalc(model.jtype{i}, 0);
            S{i} = zeros(6, 1);
        else
            [XJ, S{i}] = jcalc(model.jtype{i}, q(i));
        end

        p = model.parent(i);
        if is_flexible_link(model, p)
            edgeB{i} = endpoint_shape(model, i, p, B{p}, ...
                numel(etaDotState{p}));
            XdefParent = local_deformation_xform(edgeB{i}, etaState{p});
        else
            edgeB{i} = zeros(6, 0);
            XdefParent = eye(6);
        end
        Xup{i} = XJ * XdefParent * model.Xtree{i};
        vParent = Xup{i} * v{p};
        edgeR{i} = spatial_rotation(Xup{i});
        if is_flexible_link(model, p)
            vFlexParent = edgeR{i} * edgeB{i} * etaDotState{p};
        else
            vFlexParent = zeros(6, 1);
        end

        vJ = S{i} * qd(i) + vFlexParent;
        v{i} = vParent + vJ;
        gamma{i} = crm(v{i}) * vJ;
    end

    if is_flexible_link(model, i)
        body = local_body(model, i, true);
        [Irr{i}, Ire{i}, Mee{i}, C{i}, c{i}] = ...
            floating_frame_terms(body, v{i}, eta_i, eta_d_i);
    else
        body = local_body(model, i, false);
        Irr{i} = body.get_spatial_inertia();
        Ire{i} = zeros(6, 0);
        Mee{i} = zeros(0, 0);
        C{i} = crf(v{i}) * Irr{i} * v{i};
        c{i} = zeros(0, 1);
    end

    if ~isempty(f_ext{i})
        C{i} = C{i} - f_ext{i}(:);
        if ~isempty(B{i})
            c{i} = c{i} - B{i}.' * f_ext{i}(:);
        end
    end

    % Initialize the composite-link quantities
    Irrp{i} = Irr{i};
    Irep{i} = Ire{i};
    Ierp{i} = Ire{i}';
    Meep{i} = Mee{i};
    Cp{i} = C{i};
    cp{i} = c{i};
end

% Backward pass (Pass 2: Articulated body inertia and bias forces)
for k = NB:-1:2
    if is_flexible_link(model, k)
        Kp{k} = Irrp{k} - Irep{k} * (Meep{k} \ Ierp{k});
        bp{k} = Cp{k} - Irep{k} * (Meep{k} \ cp{k});
    else
        Kp{k} = Irrp{k};
        bp{k} = Cp{k};
    end
    if fixedJoint(k)
        H_rigid{k} = 0;
        H{k} = 0;
        Kg = Kp{k};
        Wg{k} = Kg * gamma{k} + bp{k};
    else
        H_rigid{k} = S{k}' * Irrp{k} * S{k};
        H{k} = S{k}' * Kp{k} * S{k};
        E = eye(6) - S{k} * (H{k} \ (S{k}' * Kp{k}));
        Kg = Kp{k} * E;

        Wg{k} = Kg * gamma{k} + bp{k};
        Wg{k} = Wg{k} + Kp{k} * ...
            (H{k} \ S{k} * (tau(k) - S{k}' * bp{k}));
    end
    X = Xup{k};
    j = model.parent(k);

    % Composite equations expressed in the parent link frame
    Irrp{j} = Irrp{j} + X.' * Kg * X;
    Cp{j} = Cp{j} + X.' * Wg{k};
    if is_flexible_link(model, j)
        % B maps directly to the parent endpoint twist. The endpoint
        % translation is already represented by X, so the elastic side
        % of the coupling uses a rotation-only transform.
        RB = edgeR{k} * edgeB{k};
        Irep{j} = Irep{j} + X.' * Kg * RB;
        Ierp{j} = Irep{j}.';
        Meep{j} = Meep{j} + RB.' * Kg * RB;
        cp{j} = cp{j} + RB.' * Wg{k};
    end
end

% Condense the floating base after all child branches have been assembled.
if is_flexible_link(model, 1)
    Kp{1} = Irrp{1} - Irep{1} * (Meep{1} \ Ierp{1});
    bp{1} = Cp{1} - Irep{1} * (Meep{1} \ cp{1});
else
    Kp{1} = Irrp{1};
    bp{1} = Cp{1};
end

% Floating-base solve and Pass 3 initialization
aLink0 = Kp{1} \ (baseWrench - bp{1});
base_accel = aLink0 + a_grav;
qdd = zeros(NB, 1);
eta_dd = cell(NB, 1);
aLink = cell(NB, 1);
aLink{1} = aLink0;
if is_flexible_link(model, 1)
    eta_dd{1} = -(Meep{1} \ (Ierp{1} * aLink0 + cp{1}));
else
    eta_dd{1} = zeros(0, 1);
end

% Forward pass (Pass 3: Accelerations)
for j = 2:NB
    parent = model.parent(j);
    if is_flexible_link(model, parent)
        aFlexParent = edgeR{j} * edgeB{j} * eta_dd{parent};
    else
        aFlexParent = zeros(6, 1);
    end
    aParent_j = Xup{j} * aLink{parent} + aFlexParent;

    if fixedJoint(j)
        qdd(j) = 0;
    else
        qdd(j) = H{j} \ ...
            (tau(j) - S{j}' * ...
            (bp{j} + Kp{j} * (aParent_j + gamma{j})));
    end
    n{j} = S{j} * qdd(j) + gamma{j};
    aLink{j} = aParent_j + n{j};

    if is_flexible_link(model, j)
        aRigid = aLink{j};
        eta_dd{j} = -(Meep{j} \ (Ierp{j} * aRigid + cp{j}));
    else
        eta_dd{j} = zeros(0, 1);
    end
end

end

function value = local_base_field(baseState, name, defaultValue)
if isfield(baseState, name) && ~isempty(baseState.(name))
    value = baseState.(name);
else
    value = defaultValue;
end
end

