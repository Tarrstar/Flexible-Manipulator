function [tout, Xout, info] = generalized_alpha_id_integrator(ID_d, force, tspan, x0, options)
% GENERALIZED_ALPHA_ID_INTEGRATOR
% Generalized-alpha integration based on inverse dynamics (ID).
% ID_d provides generalized forces and first-order Jacobians.

%% Stage 1: Initialization

if nargin < 5 || isempty(options)
    options = struct();
end

x0 = x0(:);
n = numel(x0)/2;
e = x0(1:n);
ed = x0(n+1:end);
tspan = tspan(:);

maxStep = get_option(options,'MaxStep',min(diff(tspan)));
absTol = get_option(options,'AbsTol',1e-10);
relTol = get_option(options,'RelTol',1e-8);
maxIterations = get_option(options,'MaxIterations',30);

if numel(tspan) == 2
    nSteps = step_count(tspan(end)-tspan(1),maxStep);
    tout = linspace(tspan(1),tspan(end),nSteps+1).';
else
    tout = tspan;
end

nSteps = numel(tout)-1;


%% Stage 2: Configure generalized-alpha parameters

rho_inf = get_option(options,'SpectralRadius',0.9);

assert(isscalar(rho_inf) && isfinite(rho_inf) && ...
    rho_inf >= 0 && rho_inf <= 1, ...
    'SpectralRadius must lie in [0,1].');

alpha_m = (2*rho_inf-1)/(rho_inf+1);
alpha_f = rho_inf/(rho_inf+1);
gamma = 0.5+alpha_f-alpha_m;
beta = 0.25*(0.5+gamma)^2;
kappa = (1-alpha_f)/(1-alpha_m);

% Jacobians are precomputed in ID_d
info.jacobianSetupSeconds = 0;


%% Stage 3: Compute initial acceleration and initialize outputs

initializationStart = tic;

[edd,ok] = newton_raphson_iter( ...
    ID_d,force(tout(1)),zeros(n,1),e,ed,0,0, ...
    maxIterations,absTol,relTol);

assert(ok, 'generalized_alpha_id_integrator:InitialBalance', ...
    'Initial acceleration did not satisfy ID-tau=0.');

% Physical and filtered accelerations
a = edd;

% Initialize outputs
Xout = zeros(numel(tout),2*n);
Xout(1,:) = x0.';

info.acceleration = zeros(numel(tout),n);
info.acceleration(1,:) = edd.';

% Initialize solver statistics
info.internalSteps = 0;
info.totalIterations = 0;
info.maximumIterations = 0;

info.absoluteResidualHistory = cell(nSteps,1);
info.relativeResidualHistory = cell(nSteps,1);

info.stepStartTimes = zeros(nSteps,1);
info.stepSizes = zeros(nSteps,1);

info.initializationSeconds = toc(initializationStart);


%% Stage 4: Perform generalized-alpha time integration

runStart = tic;

for i = 1:nSteps

    t = tout(i);
    h = tout(i+1)-t;

    % Predictor
    ed_pred = ed+h*(1-gamma)*a;
    e_pred = e+h*ed+0.5*h^2*(1-2*beta)*a;
    a_pred = (alpha_f*edd-alpha_m*a)/(1-alpha_m);

    % Generalized-alpha coefficients
    cq = beta*h^2*kappa;
    cv = gamma*h*kappa;

    qBase = e_pred+beta*h^2*a_pred;
    vBase = ed_pred+gamma*h*a_pred;

    % Newton correction
    [edd_new,ok,iter,absolute,relative] = ...
        newton_raphson_iter( ...
        ID_d,force(t+h),edd,qBase,vBase,cq,cv, ...
        maxIterations,absTol,relTol);

    if ~ok
        error('generalized_alpha_id_integrator:NoConvergence', ...
            ['ID-tau Newton failed at t=%.9g with h=%.3e; ', ...
             'residual=%.3e.'], ...
            t,h,absolute(end));
    end

    % Corrector
    a = a_pred+kappa*edd_new;
    e = e_pred+h^2*beta*a;
    ed = ed_pred+h*gamma*a;
    edd = edd_new;

    % Store results
    Xout(i+1,:) = [e;ed].';
    info.acceleration(i+1,:) = edd.';

    % Update solver statistics
    info.internalSteps = info.internalSteps+1;
    info.totalIterations = info.totalIterations+iter;
    info.maximumIterations = max(info.maximumIterations,iter);

    info.absoluteResidualHistory{i} = absolute;
    info.relativeResidualHistory{i} = relative;
    info.stepStartTimes(i) = t;
    info.stepSizes(i) = h;

end


%% Stage 5: Finalize integration statistics and solver information

info.integrationSeconds = toc(runStart);

info.nonlinearSolver = 'Newton';
info.jacobianMethod = 'CasADi automatic differentiation';
info.tangentAssembly = 'Chain rule from first-order ID Jacobians';
info.residual = 'ID-tau';
info.spectralRadius = rho_inf;

end


% =========================================================================
% Local function: Newton-Raphson nonlinear solver
% =========================================================================
function [edd_new,ok,iter,absolute,relative] = ...
    newton_raphson_iter(ID_d,Q,edd_new,qBase,vBase,cq,cv, ...
                        maxIterations,absTol,relTol)

absolute = zeros(maxIterations,1);
relative = zeros(maxIterations,1);
ok = false;

Q = full(Q(:));

for iter = 1:maxIterations

    % Generalized-alpha state reconstruction
    q = qBase+cq*edd_new;
    v = vBase+cv*edd_new;
    a = edd_new;

    % Inverse dynamics and first-order Jacobians
    [F,Jq,Jv,Ja] = ID_d(q,v,a);

    F = full(F);
    Jq = full(Jq);
    Jv = full(Jv);
    Ja = full(Ja);

    % Residual and effective Newton Jacobian
    Residual = F-Q;
    H = cq*Jq+cv*Jv+Ja;

    % Convergence evaluation
    absolute(iter) = norm(Residual,inf);
    scale = max([1,norm(Q,inf),norm(F,inf)]);
    relative(iter) = absolute(iter)/scale;

    if any(~isfinite(Residual)) || any(~isfinite(H(:)))
        break;
    end

    if absolute(iter) <= absTol+relTol*scale
        ok = true;
        break;
    end

    % Newton correction
    delta = -H\Residual;

    if any(~isfinite(delta))
        break;
    end

    edd_new = edd_new+delta;

end

absolute = absolute(1:iter);
relative = relative(1:iter);

end


% =========================================================================
% Local function: Read solver option
% =========================================================================
function value = get_option(options,name,default)

value = default;

if isfield(options,name) && ~isempty(options.(name))
    value = options.(name);
end

end


% =========================================================================
% Local function: Determine number of time steps
% =========================================================================
function count = step_count(interval,h)

ratio = interval/h;
count = max(1,ceil(ratio-1e-12*max(1,ratio)));

end