classdef FEM_Beam_GeometricStiffness < FEM_Beam
    %FEM_BEAM_GEOMETRICSTIFFNESS Planar beam coefficients for JSV Eq. (56).
    %
    % The class evaluates the two constant coefficient matrices
    %
    %   C = integral_0^L rho*A*(L-x)*Phi2_x.'*Phi2_x dx
    %   D = 1/2*integral_0^L rho*A*(L^2-x^2)*Phi2_x.'*Phi2_x dx
    %
    % used by DeltaK = omega^2*D - axialBaseAcceleration*C.  It implements
    % the first-order foreshortening correction in Liu and Hong, Journal of
    % Sound and Vibration 278 (2004), Eqs. (47), (48), (53), and (56).

    properties (SetAccess = private)
        geometric_stiffness_plane = 'xy'
        geometric_quadrature_order = 4
    end

    methods
        function obj = FEM_Beam_GeometricStiffness(params)
            obj@FEM_Beam(params);

            if isfield(params, 'geometric_stiffness_plane') && ...
                    ~isempty(params.geometric_stiffness_plane)
                plane = lower(char(params.geometric_stiffness_plane));
                if ~ismember(plane, {'xy', 'xz'})
                    error(['params.geometric_stiffness_plane must be ', ...
                        '''xy'' or ''xz''.']);
                end
                obj.geometric_stiffness_plane = plane;
            end

            if isfield(params, 'geometric_quadrature_order') && ...
                    ~isempty(params.geometric_quadrature_order)
                order = params.geometric_quadrature_order;
                if ~isscalar(order) || ~isfinite(order) || ...
                        order < 1 || order ~= round(order)
                    error(['params.geometric_quadrature_order must be a ', ...
                        'positive integer.']);
                end
                obj.geometric_quadrature_order = order;
            end
        end

        function [C, D] = get_geometric_stiffness_coefficients( ...
                obj, shapeMatrix, quadratureOrder)
        %GET_GEOMETRIC_STIFFNESS_COEFFICIENTS Compute paper matrices C and D.
        % shapeMatrix maps the selected flexible coordinates into the full
        % nodal FE coordinates used by FEM_Beam.

            if nargin < 3 || isempty(quadratureOrder)
                quadratureOrder = obj.geometric_quadrature_order;
            end
            if ~isscalar(quadratureOrder) || ~isfinite(quadratureOrder) || ...
                    quadratureOrder < 1 || ...
                    quadratureOrder ~= round(quadratureOrder)
                error('quadratureOrder must be a positive integer.');
            end
            if ~isnumeric(shapeMatrix) || ~isreal(shapeMatrix) || ...
                    size(shapeMatrix, 1) ~= obj.get_total_dof() || ...
                    any(~isfinite(shapeMatrix(:)))
                error(['shapeMatrix must be a finite real matrix with ', ...
                    'get_total_dof() rows.']);
            end

            nEta = size(shapeMatrix, 2);
            C = zeros(nEta);
            D = zeros(nEta);
            [gaussPoints, gaussWeights] = ...
                gauss_legendre_1d(quadratureOrder);

            L = obj.length;
            elementLength = L / obj.elements_num;
            rhoA = obj.density * obj.width * obj.thickness;

            for elementId = 1:obj.elements_num
                xLeft = (elementId - 1) * elementLength;
                for g = 1:numel(gaussPoints)
                    s = 0.5 * (gaussPoints(g) + 1);
                    x = xLeft + elementLength * s;
                    weight = 0.5 * elementLength * gaussWeights(g);
                    slopeShape = obj.get_transverse_slope_shape( ...
                        elementId, s) * shapeMatrix;
                    slopeProduct = slopeShape.' * slopeShape;

                    C = C + weight * rhoA * (L - x) * slopeProduct;
                    D = D + 0.5 * weight * rhoA * ...
                        (L^2 - x^2) * slopeProduct;
                end
            end

            C = 0.5 * (C + C.');
            D = 0.5 * (D + D.');
        end

        function data = get_geometric_stiffness_data( ...
                obj, shapeMatrix, quadratureOrder)
        %GET_GEOMETRIC_STIFFNESS_DATA Return fields accepted by FlexBody.

            if nargin < 3
                quadratureOrder = [];
            end
            [C, D] = obj.get_geometric_stiffness_coefficients( ...
                shapeMatrix, quadratureOrder);
            data.geometric_stiffening = true;
            data.geometric_C = C;
            data.geometric_D = D;
            data.axial_base_acceleration = 0;
            if strcmp(obj.geometric_stiffness_plane, 'xy')
                data.geometric_rotation_axis = [0; 0; 1];
            else
                data.geometric_rotation_axis = [0; 1; 0];
            end
        end
    end

    methods (Access = private)
        function slopeShape = get_transverse_slope_shape( ...
                obj, elementId, s)
            le = obj.length / obj.elements_num;
            dH1 = (-6 * s + 6 * s^2) / le;
            dH2 = 1 - 4 * s + 3 * s^2;
            dH3 = (6 * s - 6 * s^2) / le;
            dH4 = -2 * s + 3 * s^2;

            elementDofs = [obj.get_node_dofs(elementId), ...
                obj.get_node_dofs(elementId + 1)];
            slopeShape = zeros(1, obj.get_total_dof());

            if strcmp(obj.geometric_stiffness_plane, 'xy')
                localIds = [2, 6, 8, 12];
                coefficients = [dH1, dH2, dH3, dH4];
            else
                localIds = [3, 5, 9, 11];
                coefficients = [dH1, -dH2, dH3, -dH4];
            end
            slopeShape(elementDofs(localIds)) = coefficients;
        end
    end
end
