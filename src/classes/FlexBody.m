classdef FlexBody < handle

    properties (SetAccess = private)
        density
        mass

        q
        qdot

        u0
        S
        weights

        nGauss
        nEta

        S_bar
        S_t0
        M_ff

        K_ff
        D_ff
    end

    methods
        function obj = FlexBody(params)
        %FLEXBODY Construct a spatial FFRF flexible body.

            if nargin == 0
                return;
            end

            if ~isstruct(params)
                error('Input params must be a struct.');
            end

            requiredFields = {'density', 'u0', 'S', 'weights'};

            for k = 1:numel(requiredFields)
                if ~isfield(params, requiredFields{k})
                    error('Missing params.%s.', requiredFields{k});
                end
            end

            if ~isnumeric(params.density) || ...
                    ~isreal(params.density) || ...
                    ~isscalar(params.density) || ...
                    ~isfinite(params.density) || ...
                    params.density <= 0
                error('params.density must be a positive finite scalar.');
            end

            obj.density = params.density;

            obj.u0 = params.u0;
            obj.S = params.S;
            obj.weights = params.weights(:);

            if ~isnumeric(obj.u0) || ~isreal(obj.u0) || ...
                    size(obj.u0, 1) ~= 3 || ...
                    any(~isfinite(obj.u0(:)))
                error('params.u0 must be a finite 3 x nGauss matrix.');
            end

            if ~isnumeric(obj.S) || ~isreal(obj.S) || ...
                    any(~isfinite(obj.S(:)))
                error('params.S must be finite and real.');
            end

            % Allow one Gauss point: S is 3 x nEta.
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

            if numel(obj.weights) ~= obj.nGauss
                error('params.weights must contain nGauss entries.');
            end

            if any(~isfinite(obj.weights)) || any(obj.weights <= 0)
                error('params.weights must contain positive finite values.');
            end

            if isfield(params, 'Kff') && ~isempty(params.Kff)
                obj.K_ff = params.Kff;

                if ~isequal(size(obj.K_ff), [obj.nEta, obj.nEta])
                    error('params.Kff must have size nEta x nEta.');
                end

                if any(~isfinite(obj.K_ff(:)))
                    error('params.Kff cannot contain Inf or NaN.');
                end

                if norm(obj.K_ff - obj.K_ff.', 'fro') > 1e-10
                    error('params.Kff must be symmetric.');
                end
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

            obj.mass = obj.density * sum(obj.weights);

            obj.S_bar = zeros(3, obj.nEta);
            obj.S_t0 = zeros(3, 1);
            obj.M_ff = zeros(obj.nEta);

            for g = 1:obj.nGauss
                wg = obj.weights(g);

                u0g = obj.u0(:, g);
                Sg = obj.S(:, :, g);

                obj.S_t0 = obj.S_t0 ...
                    + obj.density * wg * u0g;

                obj.S_bar = obj.S_bar ...
                    + obj.density * wg * Sg;

                obj.M_ff = obj.M_ff ...
                    + obj.density * wg * (Sg.' * Sg);
            end

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
            Qv_f = zeros(obj.nEta, 1);
        
            for g = 1:obj.nGauss
                wg = obj.weights(g);
        
                Sg = obj.S(:, :, g);
        
                u = obj.u0(:, g) + Sg * eta;
        
                uDot = Sg * etaDot;
        
                temp = omega_tilde * omega_tilde * u ...
                     + 2 * omega_tilde * uDot;
        
                Qv_f = Qv_f ...
                    - obj.density * wg * (Sg.' * temp);
            end
        
            % Coordinate order: [theta; R; eta]
            Qv = [Qv_theta;
                  Qv_R;
                  Qv_f];
        end
    end

    methods (Access = private)
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

            % Eq. (5.60)
            S_t = obj.S_t0 + obj.S_bar * eta;

            I_theta_theta = zeros(3);
            I_theta_f = zeros(3, obj.nEta);
            I_theta_theta_dot = zeros(3);

            for g = 1:obj.nGauss
                wg = obj.weights(g);

                Sg = obj.S(:, :, g);
                ug = obj.u0(:, g) + Sg * eta;
                ug_dot = Sg * etaDot;

                ug_tilde = get_skew_sym_matrix(ug);
                ug_dot_tilde = get_skew_sym_matrix(ug_dot);

                % Eq. (5.68)
                I_theta_theta = I_theta_theta ...
                    + obj.density * wg * ...
                    (ug_tilde.' * ug_tilde);

                I_theta_theta_dot = I_theta_theta_dot ...
                    + obj.density * wg * ...
                    (ug_dot_tilde.' * ug_tilde ...
                     + ug_tilde.' * ug_dot_tilde);

                % Eq. (5.73)
                I_theta_f = I_theta_f ...
                    + obj.density * wg * ...
                    (ug_tilde * Sg);
            end

            I_theta_theta = 0.5 * ...
                (I_theta_theta + I_theta_theta.');

            I_theta_theta_dot = 0.5 * ...
                (I_theta_theta_dot + I_theta_theta_dot.');
        end
    end
end
