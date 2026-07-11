classdef FlexBody < handle

    properties (SetAccess = private)
        density
        mass

        q
        qdot

        u0
        S
        weights
        rhoWeights

        nGauss
        nEta

        S_bar
        S_t0
        M_ff
        shapeProducts

        K_ff
        D_ff
    end

    methods
        function obj = FlexBody(params)
        %FLEXBODY Construct a spatial FFRF flexible body.
        % params.S(:, :, g) maps flexible coordinates to displacement at point g.
        % It can contain FE, modal, or CB-reduced shape data.
        % Use params.weights with params.density, or params.massWeights = rho*dV.

            if nargin == 0
                return;
            end

            if ~isstruct(params)
                error('Input params must be a struct.');
            end

            requiredFields = {'u0', 'S'};

            for k = 1:numel(requiredFields)
                if ~isfield(params, requiredFields{k})
                    error('Missing params.%s.', requiredFields{k});
                end
            end

            obj.u0 = params.u0;
            obj.S = params.S;

            if ~isnumeric(obj.u0) || ~isreal(obj.u0) || ...
                    size(obj.u0, 1) ~= 3 || ...
                    any(~isfinite(obj.u0(:)))
                error('params.u0 must be a finite 3 x nGauss matrix.');
            end

            if ~isnumeric(obj.S) || ~isreal(obj.S) || ...
                    any(~isfinite(obj.S(:)))
                error('params.S must be finite and real.');
            end

            % Allow one integration point: S is 3 x nEta.
            % In general, S(:, :, g) already includes N(x_g), FE/modal/CB maps, etc.
            if ndims(obj.S) == 2
                obj.S = reshape(obj.S, size(obj.S, 1), ...
                    size(obj.S, 2), 1);
            end

            if ndims(obj.S) ~= 3 || size(obj.S, 1) ~= 3
                error('params.S must have size 3 x nEta x nGauss.');
            end

            obj.nGauss = size(obj.S, 3);
            obj.nEta = size(obj.S, 2);

            if size(obj.u0, 2) ~= obj.nGauss
                error(['The number of columns of params.u0 must equal ', ...
                       'the number of integration points in params.S.']);
            end

            obj.initialize_integration_weights(params);

            if isfield(params, 'Kff') && ~isempty(params.Kff)
                obj.K_ff = params.Kff;

                if ~isequal(size(obj.K_ff), [obj.nEta, obj.nEta])
                    error('params.Kff must have size nEta x nEta.');
                end

                if any(~isfinite(obj.K_ff(:)))
                    error('params.Kff cannot contain Inf or NaN.');
                end

                obj.K_ff = 0.5 * (obj.K_ff + obj.K_ff.');

                % if norm(obj.K_ff - obj.K_ff.', 'fro') > 1e-10
                %     error('params.Kff must be symmetric.');
                % end
            else
                obj.K_ff = zeros(obj.nEta);
            end

            if isfield(params, 'Dff') && ~isempty(params.Dff)
                obj.D_ff = params.Dff;

                if ~isequal(size(obj.D_ff), [obj.nEta, obj.nEta])
                    error('params.Dff must have size nEta x nEta.');
                end

                if any(~isfinite(obj.D_ff(:)))
                    error('params.Dff cannot contain Inf or NaN.');
                end

                if norm(obj.D_ff - obj.D_ff.', 'fro') > 1e-10
                    error('params.Dff must be symmetric.');
                end
            else
                obj.D_ff = zeros(obj.nEta);
            end

            obj.pre_calculate();

            nq = 7 + obj.nEta;

            if isfield(params, 'q') && ~isempty(params.q)
                q0 = params.q(:);
            else
                % q = [theta; R; eta]
                q0 = [1;
                      0;
                      0;
                      0;
                      zeros(3, 1);
                      zeros(obj.nEta, 1)];
            end

            if isfield(params, 'qdot') && ~isempty(params.qdot)
                qdot0 = params.qdot(:);
            else
                qdot0 = zeros(nq, 1);
            end

            obj.set_state(q0, qdot0);
        end

        function set_state(obj, q, qdot)
        %SET_STATE Set q = [theta; R; eta] and qdot.

            nq = 7 + obj.nEta;

            q = q(:);
            qdot = qdot(:);

            if numel(q) ~= nq
                error('q must have %d elements.', nq);
            end

            if numel(qdot) ~= nq
                error('qdot must have %d elements.', nq);
            end

            if any(~isfinite(q)) || any(~isfinite(qdot))
                error('q and qdot cannot contain Inf or NaN.');
            end

            theta = q(1:4);
            thetaDot = qdot(1:4);

            if abs(norm(theta) - 1) > 1e-8
                error(['Quaternion theta = q(1:4) must satisfy ', ...
                       'norm(theta) = 1.']);
            end

            if abs(theta.' * thetaDot) > 1e-7 * max(1, norm(thetaDot))
                warning('FlexBody:QuaternionRateConstraint', ...
                    ['theta.'' * thetaDot is not approximately zero. ', ...
                     'The quaternion-rate constraint may be violated.']);
            end

            obj.q = q;
            obj.qdot = qdot;
        end

        function pre_calculate(obj)

            obj.mass = sum(obj.rhoWeights);

            obj.S_bar = zeros(3, obj.nEta);
            obj.S_t0 = zeros(3, 1);

            nAug = obj.nEta + 1;
            obj.shapeProducts = cell(3, 3);

            for k = 1:3
                for l = 1:3
                    obj.shapeProducts{k, l} = zeros(nAug);
                end
            end

            for g = 1:obj.nGauss
                wg = obj.rhoWeights(g);

                u0g = obj.u0(:, g);
                Sg = obj.S(:, :, g);
                Hg = [u0g, Sg];

                obj.S_t0 = obj.S_t0 + wg * u0g;
                obj.S_bar = obj.S_bar + wg * Sg;

                for k = 1:3
                    for l = 1:3
                        obj.shapeProducts{k, l} = ...
                            obj.shapeProducts{k, l} ...
                            + wg * (Hg(k, :).' * Hg(l, :));
                    end
                end
            end

            obj.M_ff = ...
                obj.shapeProducts{1, 1}(2:end, 2:end) ...
              + obj.shapeProducts{2, 2}(2:end, 2:end) ...
              + obj.shapeProducts{3, 3}(2:end, 2:end);

            obj.M_ff = 0.5 * (obj.M_ff + obj.M_ff.');
        end

        function M = get_mass_matrix(obj, q)
            if nargin < 2 || isempty(q)
                q = obj.q;
            end

            [theta, ~, eta] = obj.split_q(q);

            A = get_rotation_matrix_quaternion(theta);
            [~, G_bar] = get_quaternion_G_matrices(theta);

            [S_t, I_theta_theta, I_theta_f] = ...
                obj.get_dynamic_integrals(eta);

            % Eq. (5.52)
            m_RR = obj.mass * eye(3);

            % Eq. (5.55)
            m_Rtheta = -A * get_skew_sym_matrix(S_t) * G_bar;

            % Eq. (5.64)
            m_Rf = A * obj.S_bar;

            % Eq. (5.67)
            m_theta_theta = G_bar.' * I_theta_theta * G_bar;

            % Eq. (5.72)
            m_theta_f = G_bar.' * I_theta_f;

            % Eq. (5.77)
            m_ff = obj.M_ff;

            % Coordinate order: [theta; R; eta]
            M = [m_theta_theta, m_Rtheta.',  m_theta_f;
                 m_Rtheta,      m_RR,        m_Rf;
                 m_theta_f.',   m_Rf.',      m_ff];

        end

        function K = get_stiffness_matrix(obj)

            nq = 7 + obj.nEta;

            K = zeros(nq);
            K(8:end, 8:end) = obj.K_ff;
        end

        function Qe = get_elastic_force(obj, q, qdot)

            if nargin < 2 || isempty(q)
                q = obj.q;
            end

            if nargin < 3 || isempty(qdot)
                qdot = obj.qdot;
            end

            [~, ~, eta] = obj.split_q(q);

            qdot = qdot(:);

            if numel(qdot) ~= 7 + obj.nEta
                error('qdot has an incorrect dimension.');
            end

            etaDot = qdot(8:end);

            Qe = [zeros(7, 1);
                  obj.K_ff * eta + obj.D_ff * etaDot];
        end

        function Qv = get_quadratic_velocity_vec(obj, q, qdot)

            if nargin < 2 || isempty(q)
                q = obj.q;
            end
        
            if nargin < 3 || isempty(qdot)
                qdot = obj.qdot;
            end
        
            q = q(:);
            qdot = qdot(:);
        
            nq = 7 + obj.nEta;
        
            if numel(q) ~= nq
                error('q must have %d elements.', nq);
            end
        
            if numel(qdot) ~= nq
                error('qdot must have %d elements.', nq);
            end
        
            if any(~isfinite(q)) || any(~isfinite(qdot))
                error('q and qdot cannot contain Inf or NaN.');
            end
        
            [theta, ~, eta] = obj.split_q(q);
        
            thetaDot = qdot(1:4);
            etaDot   = qdot(8:end);
        
            A = get_rotation_matrix_quaternion(theta);
        
            [~, G_bar] = get_quaternion_G_matrices(theta);
            [~, G_bar_dot] = get_quaternion_G_dot_matrices(thetaDot);
        
            omega_bar = G_bar * thetaDot;
        
            omega_tilde = get_skew_sym_matrix(omega_bar);
        
            [S_t, I_theta_theta, I_theta_f, I_theta_theta_dot] = ...
                obj.get_dynamic_integrals(eta, etaDot);
        
            % Eq. (5.149), translational part
            Qv_R = -A * ( ...
                omega_tilde * omega_tilde * S_t ...
                + 2 * omega_tilde * obj.S_bar * etaDot );
        
            % Eq. (5.149), rotational Euler-parameter part
            Qv_theta = ...
                -2 * G_bar_dot.' * I_theta_theta * omega_bar ...
                -2 * G_bar_dot.' * I_theta_f * etaDot ...
                -G_bar.' * I_theta_theta_dot * omega_bar;
        
            % Eq. (5.149), elastic part
            qAug = [1; eta];
            qAugDot = [0; etaDot];

            Qv_f = -obj.contract_shape_products( ...
                omega_tilde * omega_tilde, qAug) ...
                - 2 * obj.contract_shape_products( ...
                omega_tilde, qAugDot);
        
            % Coordinate order: [theta; R; eta]
            Qv = [Qv_theta;
                  Qv_R;
                  Qv_f];
        end
    end

    methods (Access = private)
        function initialize_integration_weights(obj, params)

            hasMassWeights = isfield(params, 'massWeights') && ...
                ~isempty(params.massWeights);
            hasDensityWeights = isfield(params, 'density') && ...
                ~isempty(params.density) && ...
                isfield(params, 'weights') && ...
                ~isempty(params.weights);

            if hasMassWeights
                obj.rhoWeights = params.massWeights(:);
                obj.weights = obj.rhoWeights;
                obj.density = ones(obj.nGauss, 1);
            elseif hasDensityWeights
                obj.weights = params.weights(:);
                density = params.density;

                if ~(isnumeric(density) && isreal(density))
                    error('params.density must be numeric and real.');
                end

                if isscalar(density)
                    obj.density = density;
                    densityVector = density * ones(obj.nGauss, 1);
                else
                    obj.density = density(:);
                    densityVector = obj.density;
                end

                if numel(densityVector) ~= obj.nGauss
                    error(['params.density must be scalar or contain ', ...
                           'nGauss entries.']);
                end

                obj.rhoWeights = densityVector .* obj.weights;
            else
                error(['Provide either params.massWeights or both ', ...
                       'params.density and params.weights.']);
            end

            if numel(obj.weights) ~= obj.nGauss
                error('params.weights must contain nGauss entries.');
            end

            if numel(obj.rhoWeights) ~= obj.nGauss
                error('params.massWeights must contain nGauss entries.');
            end

            if any(~isfinite(obj.weights)) || any(obj.weights <= 0)
                error('Integration weights must be positive and finite.');
            end

            if any(~isfinite(obj.rhoWeights)) || any(obj.rhoWeights <= 0)
                error('Mass weights rho*dV must be positive and finite.');
            end
        end

        function y = contract_shape_products(obj, C, x)

            if ~isequal(size(C), [3, 3])
                error('C must be a 3 x 3 matrix.');
            end

            x = x(:);

            if numel(x) ~= obj.nEta + 1
                error('x must have nEta + 1 entries.');
            end

            y = zeros(obj.nEta, 1);

            for k = 1:3
                for l = 1:3
                    y = y + C(k, l) * ...
                        (obj.shapeProducts{k, l}(2:end, :) * x);
                end
            end
        end

        function [theta, R, eta] = split_q(obj, q)
        %SPLIT_Q Split q = [theta; R; eta].

            q = q(:);

            nq = 7 + obj.nEta;

            if numel(q) ~= nq
                error('q must have %d elements.', nq);
            end

            if any(~isfinite(q))
                error('q cannot contain Inf or NaN.');
            end

            theta = q(1:4);
            R = q(5:7);
            eta = q(8:end);

            if abs(norm(theta) - 1) > 1e-8
                error('q(1:4) must be a unit quaternion.');
            end
        end

        function [S_t, I_theta_theta, I_theta_f, I_theta_theta_dot] = ...
                get_dynamic_integrals(obj, eta, etaDot)

            eta = eta(:);

            if numel(eta) ~= obj.nEta
                error('eta must have nEta elements.');
            end

            if nargin < 3 || isempty(etaDot)
                etaDot = zeros(obj.nEta, 1);
            else
                etaDot = etaDot(:);

                if numel(etaDot) ~= obj.nEta
                    error('etaDot must have nEta elements.');
                end
            end

            S_t = obj.S_t0 + obj.S_bar * eta;

            qAug = [1; eta];
            qAugDot = [0; etaDot];

            J11 = obj.shapeProducts{1, 1};
            J12 = obj.shapeProducts{1, 2};
            J13 = obj.shapeProducts{1, 3};
            J21 = obj.shapeProducts{2, 1};
            J22 = obj.shapeProducts{2, 2};
            J23 = obj.shapeProducts{2, 3};
            J31 = obj.shapeProducts{3, 1};
            J32 = obj.shapeProducts{3, 2};
            J33 = obj.shapeProducts{3, 3};

            qf = @(A) qAug.' * A * qAug;
            qfdot = @(A) qAugDot.' * A * qAug ...
                + qAug.' * A * qAugDot;

            I_theta_theta = [ ...
                qf(J22 + J33), -qf(J12),       -qf(J13);
               -qf(J21),        qf(J33 + J11), -qf(J23);
               -qf(J31),       -qf(J32),        qf(J11 + J22)];

            I_theta_theta_dot = [ ...
                qfdot(J22 + J33), -qfdot(J12),       -qfdot(J13);
               -qfdot(J21),        qfdot(J33 + J11), -qfdot(J23);
               -qfdot(J31),       -qfdot(J32),        qfdot(J11 + J22)];

            I_theta_f = [ ...
               -qAug.' * J32(:, 2:end) + qAug.' * J23(:, 2:end);
                qAug.' * J31(:, 2:end) - qAug.' * J13(:, 2:end);
               -qAug.' * J21(:, 2:end) + qAug.' * J12(:, 2:end)];

            I_theta_theta = 0.5 * (I_theta_theta + I_theta_theta.');
            I_theta_theta_dot = 0.5 * ...
                (I_theta_theta_dot + I_theta_theta_dot.');
        end
    end
end
