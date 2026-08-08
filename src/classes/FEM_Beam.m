classdef FEM_Beam
    properties
        length
        width
        thickness
        density
        youngs_modulus
        poisson_ratio = 0.33
        shear_modulus = []
        elements_num
        fixed_root = false
    end

    methods
        function obj = FEM_Beam(params)
            if nargin < 1
                return;
            end

            requiredFields = { ...
                'length', 'width', 'thickness', 'density', ...
                'youngs_modulus', 'elements_num'};

            for k = 1:numel(requiredFields)
                if ~isfield(params, requiredFields{k})
                    error('Missing beam.%s.', requiredFields{k});
                end
            end

            obj.length = params.length;
            obj.width = params.width;
            obj.thickness = params.thickness;
            obj.density = params.density;
            obj.youngs_modulus = params.youngs_modulus;
            obj.elements_num = params.elements_num;

            if isfield(params, 'poisson_ratio')
                obj.poisson_ratio = params.poisson_ratio;
            end

            if isfield(params, 'shear_modulus')
                obj.shear_modulus = params.shear_modulus;
            end

            if isfield(params, 'fixed_root')
                obj.fixed_root = params.fixed_root;
            end
        end

        function [M_beam, K_beam] = get_matrices_3d(obj)
            L = obj.length;
            nElement = obj.elements_num;
            l = L / nElement;

            width = obj.width;
            thickness = obj.thickness;

            rho = obj.density;
            E = obj.youngs_modulus;

            if isempty(obj.shear_modulus)
                G = E / (2 * (1 + obj.poisson_ratio));
            else
                G = obj.shear_modulus;
            end

            A = width * thickness;

            % Local coordinate:
            % x : beam axis
            % y : width direction
            % z : thickness direction
            Izz = width * thickness^3 / 12;
            Iyy = thickness * width^3 / 12;

            J = FEM_Beam.get_rectangular_torsion_constant(width, thickness);
            Ip = Iyy + Izz;

            Ke = zeros(12, 12);
            Me = zeros(12, 12);

            Ke_u = E * A / l * ...
                [ 1, -1;
                 -1,  1];

            Me_u = rho * A * l / 6 * ...
                [2, 1;
                 1, 2];

            idx_u = [1, 7];
            Ke(idx_u, idx_u) = Ke(idx_u, idx_u) + Ke_u;
            Me(idx_u, idx_u) = Me(idx_u, idx_u) + Me_u;

            Ke_tx = G * J / l * ...
                [ 1, -1;
                 -1,  1];

            Me_tx = rho * Ip * l / 6 * ...
                [2, 1;
                 1, 2];

            idx_tx = [4, 10];
            Ke(idx_tx, idx_tx) = Ke(idx_tx, idx_tx) + Ke_tx;
            Me(idx_tx, idx_tx) = Me(idx_tx, idx_tx) + Me_tx;

            Ke_v = E * Izz / l^3 * ...
                [ 12,    6*l,   -12,    6*l;
                  6*l,  4*l^2, -6*l,  2*l^2;
                 -12,   -6*l,    12,   -6*l;
                  6*l,  2*l^2, -6*l,  4*l^2];

            Me_v = rho * A * l / 420 * ...
                [ 156,    22*l,    54,   -13*l;
                  22*l,  4*l^2,  13*l,  -3*l^2;
                   54,    13*l,   156,   -22*l;
                 -13*l, -3*l^2, -22*l,   4*l^2];

            idx_v = [2, 6, 8, 12];
            Ke(idx_v, idx_v) = Ke(idx_v, idx_v) + Ke_v;
            Me(idx_v, idx_v) = Me(idx_v, idx_v) + Me_v;

            Ke_w_standard = E * Iyy / l^3 * ...
                [ 12,    6*l,   -12,    6*l;
                  6*l,  4*l^2, -6*l,  2*l^2;
                 -12,   -6*l,    12,   -6*l;
                  6*l,  2*l^2, -6*l,  4*l^2];

            Me_w_standard = rho * A * l / 420 * ...
                [ 156,    22*l,    54,   -13*l;
                  22*l,  4*l^2,  13*l,  -3*l^2;
                   54,    13*l,   156,   -22*l;
                 -13*l, -3*l^2, -22*l,   4*l^2];

            Tsign = diag([1, -1, 1, -1]);

            Ke_w = Tsign' * Ke_w_standard * Tsign;
            Me_w = Tsign' * Me_w_standard * Tsign;

            idx_w = [3, 5, 9, 11];
            Ke(idx_w, idx_w) = Ke(idx_w, idx_w) + Ke_w;
            Me(idx_w, idx_w) = Me(idx_w, idx_w) + Me_w;

            nNode = nElement + 1;
            dofPerNode = 6;
            totalDof = dofPerNode * nNode;

            I_all = [];
            J_all = [];
            K_all = [];
            M_all = [];

            for e = 1:nElement
                node1 = e;
                node2 = e + 1;

                dofs = [ ...
                    obj.get_node_dofs(node1), ...
                    obj.get_node_dofs(node2)];

                [ii, jj] = ndgrid(dofs, dofs);

                I_all = [I_all; ii(:)];
                J_all = [J_all; jj(:)];

                K_all = [K_all; Ke(:)];
                M_all = [M_all; Me(:)];
            end

            K_beam = sparse(I_all, J_all, K_all, totalDof, totalDof);
            M_beam = sparse(I_all, J_all, M_all, totalDof, totalDof);

            K_beam = 0.5 * (K_beam + K_beam');
            M_beam = 0.5 * (M_beam + M_beam');

            if obj.fixed_root
                freeDofs = obj.get_free_dofs();

                K_beam = K_beam(freeDofs, freeDofs);
                M_beam = M_beam(freeDofs, freeDofs);
            end
        end

        function [feShapeMatrix, stiffnessDofs] = ...
                build_direct_fe_shape_matrix(obj)
            totalDof = obj.get_total_dof();
            freeDofs = obj.get_free_dofs();

            activeDofs = freeDofs;
            [~, stiffnessDofs] = ismember(activeDofs, freeDofs);

            if any(stiffnessDofs == 0)
                error('Could not map active FE DOFs into the stiffness matrix.');
            end

            feShapeMatrix = zeros(totalDof, numel(activeDofs));

            for i = 1:numel(activeDofs)
                feShapeMatrix(activeDofs(i), i) = 1;
            end
        end

        function [u0, S, weights] = build_ffrf_data( ...
                obj, shapeMatrix, nx, ny, nz)

            [xi_x, wi_x] = gauss_legendre_1d(nx);
            [xi_y, wi_y] = gauss_legendre_1d(ny);
            [xi_z, wi_z] = gauss_legendre_1d(nz);

            L = obj.length;
            W = obj.width;
            T = obj.thickness;
            yNodes = W / 2 * xi_y;
            zNodes = T / 2 * xi_z;

            nElements = obj.elements_num;
            elementLength = L / nElements;
            jac = elementLength / 2 * W / 2 * T / 2;

            nFlexibleDofs = size(shapeMatrix, 2);
            nGauss = nElements * nx * ny * nz;

            u0 = zeros(3, nGauss);
            S = zeros(3, nFlexibleDofs, nGauss);
            weights = zeros(nGauss, 1);

            g = 0;

            for elementId = 1:nElements
                xLeft = (elementId - 1) * elementLength;

                for ix = 1:nx
                    x = xLeft + elementLength / 2 * (xi_x(ix) + 1);

                    for iy = 1:ny
                        y = yNodes(iy);

                        for iz = 1:nz
                            z = zNodes(iz);

                            g = g + 1;

                            u0(:, g) = [x; y; z];
                            S(:, :, g) = ...
                                obj.shape_func(x, y, z) * shapeMatrix;
                            weights(g) = wi_x(ix) * wi_y(iy) * ...
                                wi_z(iz) * jac;
                        end
                    end
                end
            end
        end

        function N = shape_func(obj, x, y, z)
            L = obj.length;
            nElements = obj.elements_num;
            le = L / nElements;

            if x >= L
                elementId = nElements;
                s = 1;
            else
                elementId = floor(x / le) + 1;
                xLeft = (elementId - 1) * le;
                s = (x - xLeft) / le;
            end

            dofs1 = obj.get_node_dofs(elementId);
            dofs2 = obj.get_node_dofs(elementId + 1);
            elementDofs = [dofs1, dofs2];

            N = zeros(3, obj.get_total_dof());
            N(:, elementDofs) = ...
                FEM_Beam.get_eb_beam_3d_shape_matrix(le, s, y, z);
        end

        function dispMatrix = evaluate_displacement( ...
                obj, shapeMatrix, x, y, z)
            dispMatrix = obj.shape_func(x, y, z) * shapeMatrix;
        end

        function dofs = get_node_dofs(~, nodeId)
            firstDof = 6 * (nodeId - 1) + 1;
            dofs = firstDof:(firstDof + 5);
        end

        function totalDof = get_total_dof(obj)
            totalDof = 6 * (obj.elements_num + 1);
        end

        function freeDofs = get_free_dofs(obj)
            totalDof = obj.get_total_dof();

            if obj.fixed_root
                freeDofs = 7:totalDof;
            else
                freeDofs = 1:totalDof;
            end
        end
    end

    methods (Static, Access = protected)
        function N = get_eb_beam_3d_shape_matrix(le, s, y, z)
            N1 = 1 - s;
            N2 = s;

            H1 = 1 - 3*s^2 + 2*s^3;
            H2 = le * (s - 2*s^2 + s^3);
            H3 = 3*s^2 - 2*s^3;
            H4 = le * (-s^2 + s^3);

            dH1 = (-6*s + 6*s^2) / le;
            dH2 = 1 - 4*s + 3*s^2;
            dH3 = (6*s - 6*s^2) / le;
            dH4 = -2*s + 3*s^2;

            N = zeros(3, 12);

            N(1, 1) = N1;
            N(1, 7) = N2;

            N(1, 2)  = N(1, 2)  - y * dH1;
            N(1, 6)  = N(1, 6)  - y * dH2;
            N(1, 8)  = N(1, 8)  - y * dH3;
            N(1, 12) = N(1, 12) - y * dH4;

            N(2, 2)  = N(2, 2)  + H1;
            N(2, 6)  = N(2, 6)  + H2;
            N(2, 8)  = N(2, 8)  + H3;
            N(2, 12) = N(2, 12) + H4;

            N(1, 3)  = N(1, 3)  - z * dH1;
            N(1, 5)  = N(1, 5)  + z * dH2;
            N(1, 9)  = N(1, 9)  - z * dH3;
            N(1, 11) = N(1, 11) + z * dH4;

            N(3, 3)  = N(3, 3)  + H1;
            N(3, 5)  = N(3, 5)  - H2;
            N(3, 9)  = N(3, 9)  + H3;
            N(3, 11) = N(3, 11) - H4;

            N(2, 4)  = N(2, 4)  - z * N1;
            N(2, 10) = N(2, 10) - z * N2;

            N(3, 4)  = N(3, 4)  + y * N1;
            N(3, 10) = N(3, 10) + y * N2;
        end

        function J = get_rectangular_torsion_constant(width, thickness)
            a = max(width, thickness);
            b = min(width, thickness);

            J = a * b^3 * ...
                (1/3 - 0.21 * (b/a) * (1 - b^4 / (12*a^4)));
        end
    end
end
