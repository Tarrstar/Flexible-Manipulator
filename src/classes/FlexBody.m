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

        geometric_stiffening
        geometric_C
        geometric_D
        geometric_rotation_axis
        axial_base_acceleration
        geometric_stiffness_data_available
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

            % MATLAB drops trailing singleton dimensions, in particular
            % for a rigid FlexBody with size 3-by-0-by-1.  The 2-D empty
            % representation is nevertheless a valid one-integration-
            % point shape matrix with nEta = 0.
            if (ndims(obj.S) ~= 2 && ndims(obj.S) ~= 3) || ...
                    size(obj.S, 1) ~= 3
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

            % Optional first-order geometric-stiffening correction from
            % Liu and Hong, JSV 278 (2004), Eqs. (47), (48), (53), (56):
            %   DeltaK = omega^2*D - axialBaseAcceleration*C.
            % It is disabled by default so existing models are unchanged.
            obj.geometric_stiffening = false;
            obj.geometric_C = zeros(obj.nEta);
            obj.geometric_D = zeros(obj.nEta);
            obj.geometric_rotation_axis = [0; 0; 1];
            obj.axial_base_acceleration = 0;
            obj.geometric_stiffness_data_available = false;

            if isfield(params, 'geometric_stiffening') && ...
                    ~isempty(params.geometric_stiffening)
                obj.geometric_stiffening = ...
                    obj.validate_logical_scalar( ...
                    params.geometric_stiffening, ...
                    'params.geometric_stiffening');
            elseif isfield(params, 'geometricStiffening') && ...
                    ~isempty(params.geometricStiffening)
                obj.geometric_stiffening = ...
                    obj.validate_logical_scalar( ...
                    params.geometricStiffening, ...
                    'params.geometricStiffening');
            end

            hasGeometricC = isfield(params, 'geometric_C') && ...
                ~isempty(params.geometric_C);
            hasGeometricD = isfield(params, 'geometric_D') && ...
                ~isempty(params.geometric_D);
            if hasGeometricC ~= hasGeometricD
                error(['params.geometric_C and params.geometric_D must ', ...
                    'be provided together.']);
            end
            if hasGeometricC
                obj.geometric_C = obj.validate_geometric_matrix( ...
                    params.geometric_C, 'params.geometric_C');
                obj.geometric_D = obj.validate_geometric_matrix( ...
                    params.geometric_D, 'params.geometric_D');
                obj.geometric_stiffness_data_available = true;
            elseif obj.geometric_stiffening
                error(['Geometric stiffening is enabled, but ', ...
                    'params.geometric_C and params.geometric_D are missing.']);
            end

            if isfield(params, 'geometric_rotation_axis') && ...
                    ~isempty(params.geometric_rotation_axis)
                axis = params.geometric_rotation_axis(:);
                if numel(axis) ~= 3 || any(~isfinite(axis)) || norm(axis) == 0
                    error(['params.geometric_rotation_axis must be a ', ...
                        'finite nonzero 3-vector.']);
                end
                obj.geometric_rotation_axis = axis / norm(axis);
            end

            if isfield(params, 'axial_base_acceleration') && ...
                    ~isempty(params.axial_base_acceleration)
                obj.set_axial_base_acceleration( ...
                    params.axial_base_acceleration);
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

            % Kgeometric = obj.get_geometric_stiffness_matrix(q, qdot);

            Qe = [zeros(7, 1);
                  (obj.K_ff) * eta ...
                  + obj.D_ff * etaDot];
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
        
            if (~is_casadi_type(q) && any(~isfinite(q))) || ...
                    (~is_casadi_type(qdot) && any(~isfinite(qdot)))
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

        function [Irr, Ire, Mee, C, c] = update_matrices(obj, v, eta, eta_d)
         
            if ~is_casadi_type(v)
                validateattributes(v, {'numeric'}, {'real', 'finite', 'vector', 'numel', 6});
            end
            if ~is_casadi_type(eta)
                validateattributes(eta, {'numeric'}, {'real', 'finite', 'numel', obj.nEta});
            end
            if ~is_casadi_type(eta_d)
                validateattributes(eta_d, {'numeric'}, {'real', 'finite', 'numel', obj.nEta});
            end
            assert(numel(v) == 6 && numel(eta) == obj.nEta && numel(eta_d) == obj.nEta, ...
                'Floating-frame input dimensions are incorrect.');
            v = v(:); eta = eta(:); eta_d = eta_d(:);
            omega = v(1:3); u = v(4:6);
            W = get_skew_sym_matrix(omega); W2 = W * W;
            omegaCrossU = W * u;
            qAug = [1; eta]; qAugDot = [0; eta_d];
        
            % Reuse each product contraction for J, Jdot, H and the elastic bias.
            % shapeProducts{k,l} = integral [u0_k; S_k']*[u0_l, S_l] dm.
            zero = 0 * v(1) + 0 * qAug(1);
            secondMoment = repmat(zero, 3, 3);
            secondMomentDot = repmat(zero, 3, 3);
            shapePosition = cell(3, 3);
            c = obj.S_bar.' * omegaCrossU;
            for k = 1:3
                for l = 1:3
                    product = obj.shapeProducts{k, l};
                    positionProduct = product * qAug;
                    velocityProduct = product * qAugDot;
                    secondMoment(k, l) = qAug.' * positionProduct;
                    secondMomentDot(k, l) = qAugDot.' * positionProduct ...
                        + qAug.' * velocityProduct;
                    shapePosition{k, l} = reshape(positionProduct(2:end), obj.nEta, 1);
                    shapeVelocity = reshape(velocityProduct(2:end), obj.nEta, 1);
                    c = c + W2(k, l) * shapePosition{k, l} ...
                        + 2 * W(k, l) * shapeVelocity;
                end
            end
            J = trace(secondMoment) * eye(3) - secondMoment;
            J = 0.5 * (J + J.');
            Jdot = trace(secondMomentDot) * eye(3) - secondMomentDot;
            Jdot = 0.5 * (Jdot + Jdot.');
            H = repmat(zero, 3, obj.nEta);
            H(1, :) = reshape(shapePosition{3, 2} - shapePosition{2, 3}, 1, obj.nEta);
            H(2, :) = reshape(shapePosition{1, 3} - shapePosition{3, 1}, 1, obj.nEta);
            H(3, :) = reshape(shapePosition{2, 1} - shapePosition{1, 2}, 1, obj.nEta);
            s = obj.S_t0 + obj.S_bar * eta;
            sTilde = get_skew_sym_matrix(s);
        
            Irr = [J, sTilde; -sTilde, obj.mass * eye(3)];
            Ire = [H; obj.S_bar];
            Mee = obj.M_ff;
            % Equivalent to IrrDot*v + crf(v)*(Irr*v + Ire*eta_d).
            % IreDot*eta_d = integral (S*eta_d) x (S*eta_d) dm = 0.
            C = [Jdot * omega + W * (J * omega + H * eta_d) ...
                     + sTilde * omegaCrossU;
                 obj.mass * omegaCrossU + W2 * s ...
                     + 2 * W * (obj.S_bar * eta_d)];
            c = c + obj.K_ff * eta + obj.D_ff * eta_d;
        end
    end

    methods (Access = private)
        function value = validate_logical_scalar(~, value, name)
            if ~(islogical(value) || isnumeric(value)) || ...
                    ~isscalar(value) || ~isreal(value) || ...
                    ~isfinite(value) || ~(value == 0 || value == 1)
                error('%s must be a logical scalar.', name);
            end
            value = logical(value);
        end

        function matrix = validate_geometric_matrix(obj, matrix, name)
            if ~isnumeric(matrix) || ~isreal(matrix) || ...
                    ~isequal(size(matrix), [obj.nEta, obj.nEta]) || ...
                    any(~isfinite(matrix(:)))
                error('%s must be a finite real nEta-by-nEta matrix.', name);
            end
            matrix = 0.5 * (matrix + matrix.');
        end

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

            if ~is_casadi_type(q) && any(~isfinite(q))
                error('q cannot contain Inf or NaN.');
            end

            theta = q(1:4);
            R = q(5:7);
            eta = q(8:end);

            if ~is_casadi_type(theta) && abs(norm(theta) - 1) > 1e-8
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
