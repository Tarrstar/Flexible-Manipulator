function [tout, Xout, info] = generalized_alpha_integrator(varargin)
%GENERALIZED_ALPHA_INTEGRATOR Implicit generalized-alpha integration.
%   [T,X] = GENERALIZED_ALPHA_INTEGRATOR(RHS,TSPAN,X0,OPTIONS,ARGS...)
%   uses an ode45-like interface. State coordinates must be ordered as
%       X = [configuration; velocity].
%   RHS must return [configurationDerivative; acceleration].
%
%   OPTIONS is a struct. The standard ODE fields RelTol, AbsTol and MaxStep
%   are supported together with:
%       PositionSize       number of configuration coordinates
%       SpectralRadius     high-frequency spectral radius, default 0.9
%       MaxIterations      nonlinear iterations per step, default 30
%       MaxSubdivisions    failed-step halvings, default 8
%       MinRelaxation      lower Aitken relaxation bound, default 0.05
%       QuaternionIndices  four quaternion indices in the configuration
%
%   When QuaternionIndices is provided, angular velocity must occupy the
%   first three velocity entries and body-frame linear velocity the next
%   three. The quaternion is updated multiplicatively and translation is
%   advanced with the midpoint orientation.
%
%   The former (TEND,STEP,MODEL,FORCE,...) signature remains available for
%   existing fixed-base tests. New code should use the ode45-like form.

if nargin < 1
    error('generalized_alpha_integrator requires input arguments.');
end

if isa(varargin{1}, 'function_handle')
    [tout, Xout, info] = local_ode_interface(varargin{:});
else
    [tout, Xout, info] = local_legacy_interface(varargin{:});
end
end

function [tout, Xout, info] = local_ode_interface( ...
        rhs, tspan, x0, options, varargin)

if nargin < 3
    error('rhs, tspan, and x0 are required.');
end
if nargin < 4 || isempty(options)
    options = struct();
end
if ~isstruct(options)
    error('options must be a struct, such as one returned by odeset.');
end
if ~isvector(tspan) || numel(tspan) < 2 || any(~isfinite(tspan)) || ...
        any(diff(tspan) <= 0)
    error('tspan must be a finite, strictly increasing vector.');
end

x0 = x0(:);
if isempty(x0) || any(~isfinite(x0))
    error('x0 must be a nonempty finite state vector.');
end

opts = local_options(options, tspan, numel(x0));
[outputTimes, nominalStep] = local_output_times(tspan(:), opts.maxStep);
opts.nominalStep = nominalStep;

derivative0 = local_rhs(rhs, outputTimes(1), x0, varargin{:});
acceleration = derivative0(opts.positionSize + 1:end);

nOutput = numel(outputTimes);
tout = outputTimes;
Xout = zeros(nOutput, numel(x0));
Xout(1, :) = x0.';

state.configuration = x0(1:opts.positionSize);
state.velocity = x0(opts.positionSize + 1:end);
state.acceleration = acceleration;

rhoInf = opts.spectralRadius;
alphaM = (2 * rhoInf - 1) / (rhoInf + 1);
alphaF = rhoInf / (rhoInf + 1);
gamma = 0.5 - alphaM + alphaF;
beta = 0.25 * (1 - alphaM + alphaF)^2;
parameters = [alphaM, alphaF, gamma, beta];

statistics.internalSteps = 0;
statistics.subdividedSteps = 0;
statistics.totalIterations = 0;
statistics.maximumIterations = 0;
statistics.absoluteResidualHistory = cell(0, 1);
statistics.relativeResidualHistory = cell(0, 1);
statistics.stepStartTimes = zeros(0, 1);
statistics.stepSizes = zeros(0, 1);

for outputIndex = 2:nOutput
    targetTime = outputTimes(outputIndex);
    t = outputTimes(outputIndex - 1);
    while t < targetTime - 10 * eps(max(1, abs(targetTime)))
        dt = min(opts.maxStep, targetTime - t);
        [state, stepStats] = local_advance( ...
            rhs, varargin, t, dt, state, 0, parameters, opts);
        t = t + dt;
        statistics = local_accumulate_statistics(statistics, stepStats);
    end
    Xout(outputIndex, :) = [state.configuration; state.velocity].';
end

info = statistics;
info.numberOfOutputSteps = nOutput - 1;
info.spectralRadius = rhoInf;
if statistics.internalSteps > 0
    info.averageIterations = statistics.totalIterations / ...
        statistics.internalSteps;
else
    info.averageIterations = 0;
end
end

function [stateOut, stats] = local_advance( ...
        rhs, rhsArgs, t, dt, stateIn, depth, parameters, opts)

[stateOut, converged, iterations, absoluteResiduals, relativeResiduals] = local_single_step( ...
    rhs, rhsArgs, t, dt, stateIn, parameters, opts);

stats.internalSteps = 1;
stats.subdividedSteps = 0;
stats.totalIterations = iterations;
stats.maximumIterations = iterations;
stats.absoluteResidualHistory = {absoluteResiduals};
stats.relativeResidualHistory = {relativeResiduals};
stats.stepStartTimes = t;
stats.stepSizes = dt;

if converged
    return;
end
if depth >= opts.maxSubdivisions
    error('generalized_alpha_integrator:NotConverged', ...
        ['Nonlinear iteration failed at t=%.9g s with dt=%.3e s ', ...
         'after %d iterations and %d subdivisions.'], ...
        t, dt, iterations, depth);
end

halfStep = 0.5 * dt;
[midState, firstStats] = local_advance( ...
    rhs, rhsArgs, t, halfStep, stateIn, depth + 1, parameters, opts);
[stateOut, secondStats] = local_advance( ...
    rhs, rhsArgs, t + halfStep, halfStep, midState, depth + 1, ...
    parameters, opts);

stats.internalSteps = firstStats.internalSteps + secondStats.internalSteps;
stats.subdividedSteps = firstStats.subdividedSteps + ...
    secondStats.subdividedSteps + 1;
stats.totalIterations = firstStats.totalIterations + ...
    secondStats.totalIterations;
stats.maximumIterations = max(firstStats.maximumIterations, ...
    secondStats.maximumIterations);
stats.absoluteResidualHistory = [firstStats.absoluteResidualHistory; ...
    secondStats.absoluteResidualHistory];
stats.relativeResidualHistory = [firstStats.relativeResidualHistory; ...
    secondStats.relativeResidualHistory];
stats.stepStartTimes = [firstStats.stepStartTimes; ...
    secondStats.stepStartTimes];
stats.stepSizes = [firstStats.stepSizes; secondStats.stepSizes];
end

function [stateOut, converged, iteration, absoluteResiduals, relativeResiduals] = local_single_step( ...
        rhs, rhsArgs, t, dt, state, parameters, opts)

alphaM = parameters(1);
alphaF = parameters(2);
gamma = parameters(3);
beta = parameters(4);

configurationN = state.configuration;
velocityN = state.velocity;
accelerationN = state.acceleration;

velocityPredictor = velocityN + dt * (1 - gamma) * accelerationN;
incrementPredictor = dt * velocityN + ...
    dt^2 * (0.5 - beta) * accelerationN;
accelerationGuess = accelerationN;
previousResidual = [];
relaxation = 1;
converged = false;
absoluteResiduals = zeros(opts.maxIterations, 1);
relativeResiduals = zeros(opts.maxIterations, 1);

for iteration = 1:opts.maxIterations
    incrementNext = incrementPredictor + beta * dt^2 * accelerationGuess;
    velocityNext = velocityPredictor + gamma * dt * accelerationGuess;
    configurationAlpha = local_update_configuration( ...
        configurationN, (1 - alphaF) * incrementNext, opts);
    velocityAlpha = alphaF * velocityN + (1 - alphaF) * velocityNext;
    tAlpha = alphaF * t + (1 - alphaF) * (t + dt);

    xAlpha = [configurationAlpha; velocityAlpha];
    derivativeAlpha = local_rhs(rhs, tAlpha, xAlpha, rhsArgs{:});
    dynamicAcceleration = derivativeAlpha(opts.positionSize + 1:end);
    target = (dynamicAcceleration - alphaM * accelerationN) / ...
        (1 - alphaM);
    residual = target - accelerationGuess;

    errorAbsolute = norm(residual, inf);
    errorRelative = errorAbsolute / max(1, norm(target, inf));
    absoluteResiduals(iteration) = errorAbsolute;
    relativeResiduals(iteration) = errorRelative;
    if errorAbsolute <= opts.absoluteTolerance || ...
            errorRelative <= opts.relativeTolerance
        accelerationGuess = target;
        converged = true;
        break;
    end

    if ~isempty(previousResidual)
        residualChange = residual - previousResidual;
        denominator = residualChange.' * residualChange;
        if denominator > eps * max(1, residual.' * residual)
            relaxation = -relaxation * ...
                (previousResidual.' * residualChange) / denominator;
            relaxation = min(1, max(opts.minRelaxation, relaxation));
        end
    end
    accelerationGuess = accelerationGuess + relaxation * residual;
    previousResidual = residual;
end

absoluteResiduals = absoluteResiduals(1:iteration);
relativeResiduals = relativeResiduals(1:iteration);

stateOut = state;
if ~converged
    return;
end

incrementNext = incrementPredictor + beta * dt^2 * accelerationGuess;
stateOut.configuration = local_update_configuration( ...
    configurationN, incrementNext, opts);
stateOut.velocity = velocityPredictor + gamma * dt * accelerationGuess;
stateOut.acceleration = accelerationGuess;
end

function configurationNext = local_update_configuration( ...
        configuration, increment, opts)

if isempty(opts.quaternionIndices)
    if numel(configuration) ~= numel(increment)
        error(['PositionSize must equal the velocity size when no ', ...
            'QuaternionIndices are supplied.']);
    end
    configurationNext = configuration + increment;
    return;
end

configurationNext = configuration;
quaternionIndices = opts.quaternionIndices;
quaternion = configuration(quaternionIndices);
quaternion = quaternion / norm(quaternion);

rotationIncrement = increment(opts.angularVelocityIndices);
deltaQuaternion = local_rotation_vector_quaternion(rotationIncrement);
quaternionNext = local_quaternion_product(quaternion, deltaQuaternion);
quaternionNext = quaternionNext / norm(quaternionNext);
if quaternionNext.' * quaternion < 0
    quaternionNext = -quaternionNext;
end
configurationNext(quaternionIndices) = quaternionNext;

midQuaternion = local_quaternion_product(quaternion, ...
    local_rotation_vector_quaternion(0.5 * rotationIncrement));
midRotation = get_rotation_matrix_quaternion(midQuaternion);
translationIncrement = increment(opts.linearVelocityIndices);
configurationNext(opts.translationIndices) = ...
    configuration(opts.translationIndices) + midRotation * translationIncrement;

configurationNext(opts.additivePositionIndices) = ...
    configuration(opts.additivePositionIndices) + ...
    increment(opts.additiveVelocityIndices);
end

function quaternion = local_rotation_vector_quaternion(rotationVector)
angle = norm(rotationVector);
if angle < 1e-10
    scale = 0.5 - angle^2 / 48;
else
    scale = sin(0.5 * angle) / angle;
end
quaternion = [cos(0.5 * angle); scale * rotationVector(:)];
end

function product = local_quaternion_product(left, right)
left = left(:);
right = right(:);
leftScalar = left(1);
rightScalar = right(1);
leftVector = left(2:4);
rightVector = right(2:4);
product = [leftScalar * rightScalar - leftVector.' * rightVector; ...
    leftScalar * rightVector + rightScalar * leftVector + ...
    cross(leftVector, rightVector)];
end

function derivative = local_rhs(rhs, t, x, varargin)
derivative = rhs(t, x, varargin{:});
derivative = derivative(:);
if numel(derivative) ~= numel(x) || any(~isfinite(derivative))
    error('RHS returned an invalid derivative at t=%.9g.', t);
end
end

function opts = local_options(options, tspan, stateSize)
opts.spectralRadius = local_option(options, ...
    {'SpectralRadius', 'spectralRadius'}, 0.9);
opts.maxIterations = local_option(options, ...
    {'MaxIterations', 'maxIterations'}, 30);
opts.relativeTolerance = local_option(options, ...
    {'RelTol', 'relativeTolerance'}, 1e-8);
opts.absoluteTolerance = local_option(options, ...
    {'AbsTol', 'absoluteTolerance'}, 1e-10);
opts.maxSubdivisions = local_option(options, ...
    {'MaxSubdivisions', 'maxSubdivisions'}, 8);
opts.minRelaxation = local_option(options, ...
    {'MinRelaxation', 'minRelaxation'}, 0.05);

defaultStep = min(diff(tspan));
opts.maxStep = local_option(options, {'MaxStep', 'StepSize'}, defaultStep);
if isempty(opts.maxStep) || ~isfinite(opts.maxStep)
    opts.maxStep = defaultStep;
end

positionSize = local_option(options, ...
    {'PositionSize', 'positionSize'}, []);
if isempty(positionSize)
    if mod(stateSize, 2) ~= 0
        error(['PositionSize is required when configuration and velocity ', ...
            'have different dimensions.']);
    end
    positionSize = stateSize / 2;
end
opts.positionSize = positionSize;
velocitySize = stateSize - positionSize;
if positionSize < 1 || velocitySize < 1
    error('PositionSize must split x0 into nonempty position and velocity blocks.');
end

opts.quaternionIndices = local_option(options, ...
    {'QuaternionIndices', 'QuaternionPositionIndices'}, []);
opts.quaternionIndices = opts.quaternionIndices(:);
if isempty(opts.quaternionIndices)
    opts.angularVelocityIndices = zeros(0, 1);
    opts.linearVelocityIndices = zeros(0, 1);
    opts.translationIndices = zeros(0, 1);
    opts.additivePositionIndices = (1:positionSize).';
    opts.additiveVelocityIndices = (1:velocitySize).';
else
    if numel(opts.quaternionIndices) ~= 4 || ...
            ~isequal(opts.quaternionIndices, (1:4).')
        error('QuaternionIndices currently must be 1:4.');
    end
    if positionSize ~= velocitySize + 1 || velocitySize < 6
        error(['Quaternion states require one more configuration coordinate ', ...
            'than velocity coordinates and at least six base velocities.']);
    end
    opts.angularVelocityIndices = (1:3).';
    opts.linearVelocityIndices = (4:6).';
    opts.translationIndices = (5:7).';
    opts.additivePositionIndices = (8:positionSize).';
    opts.additiveVelocityIndices = (7:velocitySize).';
end

validateattributes(opts.spectralRadius, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>=', 0, '<=', 1});
validateattributes(opts.maxIterations, {'numeric'}, ...
    {'scalar', 'integer', '>=', 1});
validateattributes(opts.relativeTolerance, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});
validateattributes(opts.absoluteTolerance, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});
validateattributes(opts.maxSubdivisions, {'numeric'}, ...
    {'scalar', 'integer', '>=', 0});
validateattributes(opts.minRelaxation, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>', 0, '<=', 1});
validateattributes(opts.maxStep, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive'});
validateattributes(opts.positionSize, {'numeric'}, ...
    {'scalar', 'integer', '>=', 1, '<', stateSize});
end

function value = local_option(options, names, defaultValue)
value = defaultValue;
for i = 1:numel(names)
    if isfield(options, names{i}) && ~isempty(options.(names{i}))
        value = options.(names{i});
        return;
    end
end
end

function [outputTimes, nominalStep] = local_output_times(tspan, maxStep)
if numel(tspan) > 2
    outputTimes = tspan;
    nominalStep = min(diff(tspan));
    return;
end

t0 = tspan(1);
tEnd = tspan(2);
nSteps = ceil((tEnd - t0) / maxStep);
outputTimes = linspace(t0, tEnd, nSteps + 1).';
nominalStep = (tEnd - t0) / nSteps;
end

function statistics = local_accumulate_statistics(statistics, stepStats)
statistics.internalSteps = statistics.internalSteps + stepStats.internalSteps;
statistics.subdividedSteps = statistics.subdividedSteps + ...
    stepStats.subdividedSteps;
statistics.totalIterations = statistics.totalIterations + ...
    stepStats.totalIterations;
statistics.maximumIterations = max(statistics.maximumIterations, ...
    stepStats.maximumIterations);
statistics.absoluteResidualHistory = ...
    [statistics.absoluteResidualHistory; stepStats.absoluteResidualHistory];
statistics.relativeResidualHistory = ...
    [statistics.relativeResidualHistory; stepStats.relativeResidualHistory];
statistics.stepStartTimes = ...
    [statistics.stepStartTimes; stepStats.stepStartTimes];
statistics.stepSizes = [statistics.stepSizes; stepStats.stepSizes];
end

function [tout, Xout, info] = local_legacy_interface( ...
        tEnd, stepSize, model, forceInput, x0, options)

if nargin < 4
    error('tEnd, stepSize, model, and forceInput are required.');
end
if nargin < 5
    x0 = [];
end
if nargin < 6
    options = struct();
end
if isstruct(x0) && nargin < 6
    options = x0;
    x0 = [];
end

[model, modalIndices, nEta] = local_modal_layout(model);
NB = model.NB;
fixedJoint = local_fixed_joints(model, NB);
legacyState = local_legacy_initial_state( ...
    model, x0, NB, nEta, modalIndices);
[q, qDot, eta, etaDot] = local_split_legacy_state(legacyState, NB, nEta);
newState = [q; eta; qDot; etaDot];

newOptions = options;
newOptions.PositionSize = NB + nEta;
newOptions.MaxStep = stepSize;
rhs = @(t, x) local_legacy_rhs( ...
    t, x, model, forceInput, NB, nEta, modalIndices, fixedJoint);
[tout, newOutput, info] = local_ode_interface( ...
    rhs, [0, tEnd], newState, newOptions);

qOutput = newOutput(:, 1:NB);
etaOutput = newOutput(:, NB + (1:nEta));
velocityOffset = NB + nEta;
qDotOutput = newOutput(:, velocityOffset + (1:NB));
etaDotOutput = newOutput(:, velocityOffset + NB + (1:nEta));
Xout = [qOutput, qDotOutput, etaOutput, etaDotOutput];
end

function dx = local_legacy_rhs( ...
        t, x, model, forceInput, NB, nEta, modalIndices, fixedJoint)

q = x(1:NB);
eta = x(NB + (1:nEta));
velocityOffset = NB + nEta;
qDot = x(velocityOffset + (1:NB));
etaDot = x(velocityOffset + NB + (1:nEta));
tau = local_force(forceInput, t, NB, fixedJoint);
[qDDot, etaDDotCell] = FD_flex( ...
    model, q, qDot, tau, eta, etaDot);
etaDDot = local_modal_vector(etaDDotCell, modalIndices, nEta);
dx = [qDot; etaDot; qDDot(:); etaDDot];
end

function [model, modalIndices, nEta] = local_modal_layout(model)
NB = model.NB;
modalIndices = cell(NB, 1);
nextIndex = 1;
hasProvidedLayout = isfield(model, 'modal_idx') && ...
    numel(model.modal_idx) >= NB;

for i = 1:NB
    if is_flexible_link(model, i)
        body = local_body(model, i, true);
        count = body.nEta;
        if hasProvidedLayout && ~isempty(model.modal_idx{i})
            indices = model.modal_idx{i}(:);
        else
            indices = (nextIndex:(nextIndex + count - 1)).';
        end
        if numel(indices) ~= count || any(indices < 1) || ...
                any(indices ~= round(indices))
            error('model.modal_idx{%d} must contain %d valid indices.', ...
                i, count);
        end
        modalIndices{i} = indices;
        if ~isempty(indices)
            nextIndex = max(nextIndex, max(indices) + 1);
        end
    else
        modalIndices{i} = zeros(0, 1);
    end
end

allIndices = vertcat(modalIndices{:});
if isempty(allIndices)
    nEta = 0;
elseif numel(unique(allIndices)) ~= numel(allIndices)
    error('model.modal_idx entries for flexible links must not overlap.');
else
    nEta = max(allIndices);
    if ~isequal(sort(allIndices), (1:nEta).')
        error('model.modal_idx must cover the contiguous range 1:nEta.');
    end
end
model.modal_idx = modalIndices;
end

function vector = local_modal_vector(values, modalIndices, nEta)
vector = zeros(nEta, 1);
if nEta == 0
    return;
end
if ~iscell(values)
    vector = values(:);
    if numel(vector) ~= nEta
        error('FD_flex returned an invalid modal acceleration vector.');
    end
    return;
end
for i = 1:numel(modalIndices)
    indices = modalIndices{i};
    if ~isempty(indices)
        vector(indices) = values{i}(:);
    end
end
end

function tau = local_force(forceInput, t, NB, fixedJoint)
if isa(forceInput, 'function_handle')
    value = forceInput(t);
else
    value = forceInput;
end
if isscalar(value)
    tau = zeros(NB, 1);
    tau(~fixedJoint) = value;
else
    tau = value(:);
end
if numel(tau) ~= NB || any(~isfinite(tau))
    error('forceInput returned an invalid torque vector at t=%.9g.', t);
end
end

function x = local_legacy_initial_state(model, x0, NB, nEta, modalIndices)
stateSize = 2 * NB + 2 * nEta;
if ~isempty(x0)
    x = x0(:);
    if numel(x) ~= stateSize
        error('x0 must contain %d entries [q;qDot;eta;etaDot].', stateSize);
    end
elseif all(isfield(model, {'q', 'qd', 'eta', 'eta_d'}))
    eta = local_state_modal_vector(model.eta, modalIndices, nEta);
    etaDot = local_state_modal_vector(model.eta_d, modalIndices, nEta);
    x = [model.q(:); model.qd(:); eta; etaDot];
else
    x = zeros(stateSize, 1);
end
end

function vector = local_state_modal_vector(values, modalIndices, nEta)
if iscell(values)
    vector = zeros(nEta, 1);
    for i = 1:numel(modalIndices)
        if ~isempty(modalIndices{i})
            vector(modalIndices{i}) = values{i}(:);
        end
    end
else
    vector = values(:);
end
end

function [q, qDot, eta, etaDot] = local_split_legacy_state(x, NB, nEta)
q = x(1:NB);
qDot = x(NB + (1:NB));
eta = x(2 * NB + (1:nEta));
etaDot = x(2 * NB + nEta + (1:nEta));
end

function fixedJoint = local_fixed_joints(model, NB)
fixedJoint = false(NB, 1);
if ~isfield(model, 'fixed_joints') || isempty(model.fixed_joints)
    return;
end
specification = model.fixed_joints;
if islogical(specification) && numel(specification) == NB
    fixedJoint = specification(:);
else
    indices = specification(:);
    if any(indices < 1 | indices > NB | indices ~= round(indices))
        error('model.fixed_joints contains invalid indices.');
    end
    fixedJoint(indices) = true;
end
end
