function [qdd, eta_dd] = FD_flex(model, q, qd, tau, eta, eta_d, f_ext)

NB = model.NB;
if ~isfield(model, 'link_type')
    error('model.link_type must specify ''rigid'' or ''flexible'' for every link.');
end
q = q(:); qd = qd(:); tau = tau(:);
if numel(q) ~= NB || numel(qd) ~= NB || numel(tau) ~= NB
    error('q, qd and tau must contain one entry per joint.');
end

fixedJoint = false(NB, 1);
if isfield(model, 'fixed_joints') && ~isempty(model.fixed_joints)
    fixedSpec = model.fixed_joints;
    if islogical(fixedSpec) && numel(fixedSpec) == NB
        fixedJoint = fixedSpec(:);
    else
        fixedIds = fixedSpec(:);
        if any(fixedIds < 1 | fixedIds > NB | fixedIds ~= round(fixedIds))
            error('model.fixed_joints must contain valid joint indices.');
        end
        fixedJoint(fixedIds) = true;
    end
end
if nargin < 7 || isempty(f_ext)
    f_ext = cell(NB, 1);
elseif ~iscell(f_ext)
    if NB ~= 1
        error('For multiple links, f_ext must be a cell array with one force per link.');
    end
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
elseif numel(a_grav) ~= 6
    error('model.gravity must contain either 3 linear or 6 spatial components.');
end

% Forward pass (Pass 1: Kinematics and local dynamics terms)
Xup = cell(NB, 1); S = cell(NB, 1); B = cell(NB, 1);
etaState = cell(NB, 1); etaDotState = cell(NB, 1);
edgeB = cell(NB, 1); edgeR = cell(NB, 1);
v = cell(NB, 1); n = cell(NB, 1); gamma = cell(NB, 1);
Irr = cell(NB, 1); Ire = cell(NB, 1); Mee = cell(NB, 1);
C = cell(NB, 1); c = cell(NB, 1);
Irrp = cell(NB, 1); Irep = cell(NB, 1); Ierp = cell(NB, 1); Meep = cell(NB, 1);
Cp = cell(NB, 1); cp = cell(NB, 1);
Kp = cell(NB, 1); bp = cell(NB, 1);
Wg = cell(NB, 1); H = cell(NB, 1); H_rigid = cell(NB, 1);

for i = 1:NB
    % 1. 获取当前连杆 i 本身的模态状态、振型矩阵 B{i} 及变形位姿 Xdef{i}
    if is_flexible_link(model, i)
        [B{i}, eta_i, eta_d_i] = local_modal_state(model, i, eta, eta_d);
    else
        B{i} = zeros(6, 0);
        eta_i = zeros(0, 1); eta_d_i = zeros(0, 1);
    end
    etaState{i} = eta_i;
    etaDotState{i} = eta_d_i;

    % 2. 铰接运动子空间及包含模态耦合的相对速度 vJ
    if fixedJoint(i)
        [XJ, ~] = jcalc(model.jtype{i}, 0);
        S{i} = zeros(6, 1);
    else
        [XJ, S{i}] = jcalc(model.jtype{i}, q(i));
    end

    % 3. 父连杆 (Parent) 几何变形传递与 Xup 计算
    p = model.parent(i);
    if p ~= 0
        if is_flexible_link(model, p)
            edgeB{i} = endpoint_shape(model, i, p, B{p}, ...
                numel(etaDotState{p}));
            XdefParent = local_deformation_xform(edgeB{i}, etaState{p});
        else
            edgeB{i} = zeros(6, 0);
            XdefParent = eye(6);
        end
        Xup{i} = XJ * XdefParent * model.Xtree{i};
        % 跨连杆变换矩阵：铰接变换 * 父连杆末端变形 * 树拓扑偏移
        vParent = Xup{i} * v{p};
        edgeR{i} = spatial_rotation(Xup{i});
        if is_flexible_link(model, p)
            vFlexParent = edgeR{i} * edgeB{i} * etaDotState{p};
        else
            vFlexParent = zeros(6, 1);
        end
    else
        % 基座连杆 (Base)
        Xup{i} = XJ * model.Xtree{i};
        vParent = zeros(6, 1);
        edgeR{i} = spatial_rotation(Xup{i});
        edgeB{i} = zeros(6, 0);
        vFlexParent = zeros(6, 1);
    end

    % 4. 当前连杆 i 的总空间速度与偏置加速度 gamma
    vJ = S{i} * qd(i) + vFlexParent;
    v{i} = vParent + vJ;
    if ~isequal(size(v{i}), [6, 1])
        error(['Link %d velocity must be 6-by-1; vParent is %d-by-%d, ', ...
            'joint/endpoint velocity is %d-by-%d.'], i, ...
            size(vParent, 1), size(vParent, 2), size(vJ, 1), size(vJ, 2));
    end
    gamma{i} = crm(v{i}) * vJ;

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
for k = NB:-1:1
    if is_flexible_link(model, k)
        if isempty(Meep{k})
            error('Flexible link %d has no composite modal mass matrix.', k);
        end
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

    if model.parent(k) == 0
        continue
    end
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

% Forward pass (Pass 3: Accelerations)
baseAccel = -a_grav;
qdd = zeros(NB, 1);
eta_dd = cell(NB, 1);
aLink = cell(NB, 1);

for j = 1:NB
    parent = model.parent(j);
    if parent == 0
        aParent_j = Xup{j} * baseAccel;
    else
        % 父连杆的总空间加速度 = 父连杆参考系加速度 + 父连杆模态加速度
        
        % 利用 Xup{j} 将父节点总加速度映射到当前连杆 j 坐标系
        aParent_j = Xup{j} * aLink{parent};
        if is_flexible_link(model, parent)
            aFlexParent = edgeR{j} * edgeB{j} * eta_dd{parent};
        else
            aFlexParent = zeros(6, 1);
        end
        aParent_j = Xup{j} * aLink{parent} + aFlexParent;
    end
    
    % 求解关节加速度 qdd
    if fixedJoint(j)
        qdd(j) = 0;
    else
        qdd(j) = H{j} \ ...
            (tau(j) - S{j}' * ...
            (bp{j} + Kp{j} * (aParent_j + gamma{j})));
    end
    n{j} = S{j} * qdd(j) + gamma{j};
    aLink{j} = aParent_j + n{j};
    
    % 求解模态加速度 eta_dd
    if is_flexible_link(model, j)
        aRigid = aLink{j};
        eta_dd{j} = -(Meep{j} \ (Ierp{j} * aRigid + cp{j}));
    else
        eta_dd{j} = zeros(0, 1);
    end
end

end

function R6 = spatial_rotation(X)
R = X(1:3, 1:3);
R6 = blkdiag(R, R);
end

function Phi = endpoint_shape(model, child, parent, defaultPhi, nEta)
% Use an edge-specific endpoint map for branched flexible structures.
% Existing serial models continue to use model.Phi{parent}.
Phi = defaultPhi;
if isfield(model, 'Phi_to_child') && numel(model.Phi_to_child) >= child && ...
        ~isempty(model.Phi_to_child{child})
    Phi = model.Phi_to_child{child};
end

if ~is_flexible_link(model, parent)
    Phi = zeros(6, 0);
    return;
end
if size(Phi, 1) ~= 6 || size(Phi, 2) ~= nEta
    error(['Endpoint shape map for parent link %d and child link %d must ', ...
        'have size 6-by-%d.'], parent, child, nEta);
end
end
