function [tau, eta_dd] = ID_flex(model, q, qd, qdd, eta, eta_d, f_ext)
%ID_FLEX Inverse dynamics for a tree of rigid/flexible links.
%
% This is the inverse-dynamics specialization of Sec. 6.2 in Khalil et
% al. (2017): all active joint accelerations are prescribed, while the
% elastic accelerations are calculated.  The recursion is, in the order
% used in the paper,
%
%   1) forward initialization of kinematics and link matrices,
%   2) backward construction of composite K+, b+ and M+ee,
%   3) forward calculation of elastic accelerations, link forces and
%      active-joint torques.
%
% A flexible link must provide a FlexBody object in model.flexbody{i}; a
% rigid link must provide a RigidBody object in model.rigidbody{i}.  Phi{i}
% is the spatial elastic velocity/acceleration map (6-by-ne) at the link
% connection.  Inertial and nonlinear terms are obtained from the
% corresponding body class.

NB = model.NB;
if ~isfield(model, 'link_type')
    error('model.link_type must specify ''rigid'' or ''flexible'' for every link.');
end
q = q(:); qd = qd(:); qdd = qdd(:);
if numel(q) ~= NB || numel(qd) ~= NB || numel(qdd) ~= NB
    error('q, qd and qdd must contain one entry per joint.');
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
% elseif exist('get_gravity', 'file') == 2
%     a_grav = get_gravity(model);
else
    a_grav = zeros(6, 1);
end
if numel(a_grav) == 3
    a_grav = [0; 0; 0; a_grav];
elseif numel(a_grav) ~= 6
    error('model.gravity must contain either 3 linear or 6 spatial components.');
end

% Forward pass (paper 6.2(i)).
Xup = cell(NB, 1); S = cell(NB, 1); B = cell(NB, 1);
v = cell(NB, 1); g = cell(NB, 1); n = cell(NB, 1);
Irr = cell(NB, 1); Ire = cell(NB, 1); Ier = cell(NB, 1); Mee = cell(NB, 1);
C = cell(NB, 1); c = cell(NB, 1);
Irrp = cell(NB, 1); Irep = cell(NB, 1); Ierp = cell(NB, 1); Meep = cell(NB, 1);
Cp = cell(NB, 1); cp = cell(NB, 1);
Kp = cell(NB, 1); bp = cell(NB, 1);
Wg = cell(NB, 1);

for i = 1:NB
    [XJ, S{i}] = jcalc(model.jtype{i}, q(i));
    [B{i}, eta_i, eta_d_i] = local_modal_state(model, i, eta, eta_d);
    Xup{i} = local_deformation_xform(B{i}, eta_i) * XJ * model.Xtree{i};

    vJ = S{i} * qd(i);
    vF = B{i} * eta_d_i;
    vParent = zeros(6, 1);
    if model.parent(i) ~= 0
        vParent = Xup{i} * v{model.parent(i)};
    end
    v{i} = vParent + vJ + vF;
    n{i} = S{i} * qdd(i) + crm(v{i}) * (vJ + vF);

    if is_flexible_link(model, i)
        body = local_body(model, i, true);
        [Irr{i}, Ire{i}, Mee{i}, C{i}, c{i}] = ...
            floating_frame_terms(body, eta_i, eta_d_i);
    else
        body = local_body(model, i, false);
        Irr{i} = body.get_spatial_inertia();
        Ire{i} = zeros(6, 0);
        Mee{i} = zeros(0, 0);
        % A RigidBody has no independent floating-frame velocity.
        C{i} = crf(v{i}) * Irr{i} * v{i};
        c{i} = zeros(0, 1);
    end
    if ~isempty(f_ext{i})
        C{i} = C{i} - f_ext{i}(:);
        if ~isempty(B{i})
            c{i} = c{i} - B{i}.' * f_ext{i}(:);
        end
    end

    % Initialize the composite-link quantities.
    Irrp{i} = Irr{i};
    Irep{i} = Ire{i};
    Ierp{i} = Ire{i}';
    Meep{i} = Mee{i};
    Cp{i} = C{i};
    cp{i} = c{i};
end

% Backward pass (paper 6.2(ii)).  Every prescribed joint acceleration uses
% Kg=K+ and Wg=K+*n+b+.
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
    Kg = Kp{k};
    Wg{k} = Kg * n{k} + bp{k};
    X = Xup{k};

    if model.parent(k) == 0
        continue
    end
    j = model.parent(k);

    % Composite equations expressed in the parent link frame.
    Irrp{j} = Irrp{j} + X.' * Kg * X;
    Cp{j} = Cp{j} + X.' * Wg{k};
    if is_flexible_link(model, j)
        XB = X * B{j};
        Irep{j} = Irep{j} + X.' * Kg * XB;
        Ierp{j} = Irep{j}.';
        Meep{j} = Meep{j} + XB.' * Kg * XB;
        cp{j} = cp{j} + XB.' * Wg{k};
    end
end

% Fixed-base inverse dynamics.  If a floating base is needed, its passive
% acceleration must be obtained from Kp{1}\bp{1} before this pass.
baseAccel = -a_grav;
tau = zeros(NB, 1);
eta_dd = cell(NB, 1);
aLink = cell(NB, 1);
F = cell(NB, 1);
for j = 1:NB
    if model.parent(j) == 0
        aParent = baseAccel;
    else
        aParent = aLink{model.parent(j)};
    end
    if is_flexible_link(model, j)
        aRigid = Xup{j} * aParent + n{j};
        eta_dd{j} = -(Meep{j} \ (Ierp{j} * aRigid + cp{j}));
    else
        eta_dd{j} = zeros(0, 1);
    end
    aLink{j} = Xup{j} * aParent + n{j} + B{j} * eta_dd{j};

    % Kp and Wg are already condensed with respect to the link modal
    % acceleration; adding B*eta_dd here would count that elimination twice.
    F{j} = Kp{j} * (Xup{j} * aParent) + Wg{j};
    tau(j) = S{j}.' * F{j};
end
end

function [Phi, eta_i, eta_d_i] = local_modal_state(model, i, eta, eta_d)
Phi = zeros(6, 0);
if isfield(model, 'Phi') && numel(model.Phi) >= i && ...
        ~isempty(model.Phi{i}) && is_flexible_link(model, i)
    Phi = model.Phi{i};
else
    body = [];
    if is_flexible_link(model, i)
        body = local_body(model, i, true);
    end
    if ~isempty(body) && body.nEta > 0
        % A zero map is appropriate for a floating frame whose modal
        % deformation does not alter the frame pose.  A nonzero endpoint
        % map can still be supplied through model.Phi{i} when a child is
        % attached away from the floating frame.
        Phi = zeros(6, body.nEta);
    end
end

n = size(Phi, 2);
if n == 0
    eta_i = zeros(0, 1); eta_d_i = zeros(0, 1); return
end
eta_i = local_link_vector(model, i, eta, n);
eta_d_i = local_link_vector(model, i, eta_d, n);
end

function tf = is_flexible_link(model, i)
types = model.link_type;
if iscell(types)
    type = types{i};
elseif isstring(types)
    if isscalar(types)
        type = types;
    else
        type = types(i);
    end
elseif ischar(types)
    type = types;
else
    error('model.link_type must be a character or string array/cell array.');
end
tf = strcmpi(strtrim(char(type)), 'flexible');
if ~tf && ~strcmpi(strtrim(char(type)), 'rigid')
    error('model.link_type{%d} must be ''rigid'' or ''flexible''.', i);
end
end

function x = local_link_vector(model, i, values, n)
if iscell(values)
    x = values{i}(:);
elseif isfield(model, 'modal_idx') && numel(model.modal_idx) >= i
    x = values(model.modal_idx{i}); x = x(:);
else
    x = values(:);
end
if numel(x) ~= n
    error('Modal coordinate length does not match size(model.Phi{%d},2).', i);
end
end

function body = local_body(model, i, isFlexible)
body = [];
if isFlexible
    names = {'flexbody', 'flex_body', 'FlexBody'};
    className = 'FlexBody';
else
    names = {'rigidbody', 'rigid_body', 'RigidBody'};
    className = 'RigidBody';
end
for name = names
    name = name{1};
    if ~isfield(model, name) || numel(model.(name)) < i
        continue
    end
    if iscell(model.(name))
        candidate = model.(name){i};
    else
        candidate = model.(name)(i);
    end
    if ~isempty(candidate)
        if isa(candidate, className)
            body = candidate;
            return
        end
    end
end
if isempty(body)
    error('Missing %s object for link %d.', className, i);
end
end

function [Irr, Ire, Mee, C, c] = floating_frame_terms(body, eta, eta_d)
% Convert FlexBody's [quaternion; translation; eta] equations to spatial
% floating-frame blocks.  The conversion uses the exact mass matrix and
% quadratic-velocity vector generated by the floating-frame formulation.
qBody = body.q;
qdotBody = body.qdot;
if numel(qBody) ~= 7 + body.nEta || numel(qdotBody) ~= 7 + body.nEta
    error('FlexBody state has an inconsistent dimension.');
end
qBody(8:end) = eta;
qdotBody(8:end) = eta_d;

M = body.get_mass_matrix(qBody);
Qv = body.get_quadratic_velocity_vec(qBody, qdotBody);
Qe = body.get_elastic_force(qBody, qdotBody);
[~, Gbar] = get_quaternion_G_matrices(qBody(1:4));
A = get_rotation_matrix_quaternion(qBody(1:4));

% omega = Gbar*theta_dot, hence theta_dot = Gbar'*omega on the
% unit-quaternion tangent space.  FlexBody stores R_dot in inertial
% coordinates, while spatial vectors use the floating-frame coordinates:
% R_dot = A*v.
Jr = [Gbar.', zeros(4, 3); zeros(3, 3), A];
Mrr = M(1:7, 1:7);
Irr = Jr.' * Mrr * Jr;
nEta = body.nEta;
if nEta == 0
    Ire = zeros(6, 0);
    Mee = zeros(0, 0);
else
    % Keep the block orientations explicit; MATLAB can collapse a
    % one-mode row/column when slicing a 7-by-(1) block.
    Ire = reshape(Jr.' * M(1:7, 8:(7 + nEta)), 6, nEta);
    Mee = reshape(M(8:(7 + nEta), 8:(7 + nEta)), nEta, nEta);
end
C = Jr.' * Qv(1:7);
c = Qv(8:end) + Qe(8:end);
end

function X = local_deformation_xform(Phi, eta)
if isempty(Phi)
    X = eye(6); return
end
[R, ~, delta] = get_deformation_transform(Phi, eta);
u = delta(4:6);
% if exist('plux', 'file') == 2
%     X = plux(R, u.');
% else
    X = [R, zeros(3); -R * get_skew_sym_matrix(u), R];
% end
end
