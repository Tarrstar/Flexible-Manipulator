classdef FEM_Beam_CB_Reduction < FEM_Beam
    properties
        reductionMatrix = []
        fullShapeMatrix = []
        fullStiffnessDofs = []
        boundaryNodeIds = []
        nInteriorModes = 0
        cbInfo = struct()
    end

    methods
        function obj = FEM_Beam_CB_Reduction(params)
            obj@FEM_Beam(params);

            if isfield(params, 'reductionMatrix')
                obj.reductionMatrix = params.reductionMatrix;
            end

            if isfield(params, 'fullShapeMatrix')
                obj.fullShapeMatrix = params.fullShapeMatrix;
            end

            if isfield(params, 'fullStiffnessDofs')
                obj.fullStiffnessDofs = params.fullStiffnessDofs;
            end

            if isfield(params, 'boundaryNodeIds')
                obj.boundaryNodeIds = params.boundaryNodeIds(:).';
            end

            if isfield(params, 'nInteriorModes')
                obj.nInteriorModes = params.nInteriorModes;
            end
        end

        function [reductionMatrix, info] = build_cb_reduction( ...
                obj, nInteriorModes, boundaryNodeIds)

            if nargin < 2 || isempty(nInteriorModes)
                nInteriorModes = obj.nInteriorModes;
            end

            if nargin < 3 || isempty(boundaryNodeIds)
                boundaryNodeIds = obj.boundaryNodeIds;
            end

            if isempty(boundaryNodeIds)
                boundaryNodeIds = [1, obj.elements_num + 1];
            end

            [MFull, KFull] = obj.get_matrices_3d();
            [~, stiffnessDofs] = obj.build_direct_fe_shape_matrix();

            MFull = full(MFull(stiffnessDofs, stiffnessDofs));
            KFull = full(KFull(stiffnessDofs, stiffnessDofs));
            MFull = 0.5 * (MFull + MFull.');
            KFull = 0.5 * (KFull + KFull.');

            freeDofs = obj.get_free_dofs();
            boundaryGlobalDofs = [];

            for k = 1:numel(boundaryNodeIds)
                boundaryGlobalDofs = [ ...
                    boundaryGlobalDofs, ...
                    obj.get_node_dofs(boundaryNodeIds(k))];
            end

            boundaryGlobalDofs = unique(boundaryGlobalDofs, 'stable');
            [isBoundaryFree, boundaryCoordIds] = ismember( ...
                boundaryGlobalDofs, freeDofs);
            boundaryCoordIds = boundaryCoordIds(isBoundaryFree);
            boundaryCoordIds = unique(boundaryCoordIds(:).', 'stable');

            nFull = numel(freeDofs);

            if isempty(boundaryCoordIds)
                error('No boundary DOFs remain in the flexible coordinate set.');
            end

            allCoordIds = 1:nFull;
            interiorCoordIds = setdiff(allCoordIds, boundaryCoordIds, ...
                'stable');

            nb = numel(boundaryCoordIds);
            ni = numel(interiorCoordIds);
            nModes = min(max(0, nInteriorModes), ni);

            reductionMatrix = zeros(nFull, nb + nModes);
            reductionMatrix(boundaryCoordIds, 1:nb) = eye(nb);

            modalFrequencies = zeros(nModes, 1);

            if ni > 0
                Kii = KFull(interiorCoordIds, interiorCoordIds);
                Kib = KFull(interiorCoordIds, boundaryCoordIds);
                Mii = MFull(interiorCoordIds, interiorCoordIds);

                constraintModes = -Kii \ Kib;
                reductionMatrix(interiorCoordIds, 1:nb) = constraintModes;

                if nModes > 0
                    [V, D] = eig(Kii, Mii);
                    lambda = real(diag(D));
                    valid = isfinite(lambda) & lambda > 0;

                    V = real(V(:, valid));
                    lambda = lambda(valid);

                    [lambda, order] = sort(lambda, 'ascend');
                    V = V(:, order);

                    if numel(lambda) < nModes
                        nModes = numel(lambda);
                        reductionMatrix = reductionMatrix(:, 1:(nb + nModes));
                    end

                    V = V(:, 1:nModes);
                    lambda = lambda(1:nModes);

                    for k = 1:nModes
                        modalMass = V(:, k).' * Mii * V(:, k);

                        if modalMass <= 0 || ~isfinite(modalMass)
                            error('Invalid fixed-interface modal mass.');
                        end

                        V(:, k) = V(:, k) / sqrt(modalMass);
                    end

                    reductionMatrix(interiorCoordIds, nb + (1:nModes)) = V;
                    modalFrequencies = sqrt(lambda) / (2 * pi);
                end
            end

            info = struct();
            info.boundaryNodeIds = boundaryNodeIds;
            info.boundaryGlobalDofs = boundaryGlobalDofs;
            info.boundaryCoordIds = boundaryCoordIds;
            info.interiorCoordIds = interiorCoordIds;
            info.nBoundaryDofs = nb;
            info.nInteriorDofs = ni;
            info.nInteriorModes = nModes;
            info.modalFrequencies = modalFrequencies;
        end

        function [shapeMatrix, stiffnessMatrix, reductionMatrix, info] = ...
                build_reduced_shape_matrix(obj, reductionSpec)

            if nargin < 2 || isempty(reductionSpec)
                reductionSpec = obj.reductionMatrix;
            end

            if isempty(reductionSpec)
                [reductionMatrix, info] = obj.build_cb_reduction();
            elseif isscalar(reductionSpec)
                [reductionMatrix, info] = obj.build_cb_reduction( ...
                    reductionSpec);
            else
                reductionMatrix = reductionSpec;
                info = struct();
            end

            if isempty(obj.fullShapeMatrix)
                [fullShapeMatrix, stiffnessDofs] = ...
                    obj.build_direct_fe_shape_matrix();
            else
                fullShapeMatrix = obj.fullShapeMatrix;
                stiffnessDofs = obj.fullStiffnessDofs;

                if isempty(stiffnessDofs)
                    stiffnessDofs = 1:size(fullShapeMatrix, 2);
                end
            end

            [~, KFull] = obj.get_matrices_3d();
            KFull = full(KFull(stiffnessDofs, stiffnessDofs));
            KFull = 0.5 * (KFull + KFull.');

            if size(reductionMatrix, 1) ~= size(fullShapeMatrix, 2)
                error(['reductionMatrix row count must match the ', ...
                       'full FE coordinate count.']);
            end

            shapeMatrix = fullShapeMatrix * reductionMatrix;
            stiffnessMatrix = reductionMatrix.' * KFull * reductionMatrix;
            stiffnessMatrix = 0.5 * ...
                (stiffnessMatrix + stiffnessMatrix.');
        end

        function [u0, S, weights, shapeMatrix, stiffnessMatrix, ...
                reductionMatrix, info] = build_ffrf_data( ...
                obj, reductionSpec, nx, ny, nz)

            if nargin < 2 || isempty(reductionSpec)
                reductionSpec = obj.reductionMatrix;
            end

            [shapeMatrix, stiffnessMatrix, reductionMatrix, info] = ...
                obj.build_reduced_shape_matrix(reductionSpec);

            [u0, S, weights] = build_ffrf_data@FEM_Beam( ...
                obj, shapeMatrix, nx, ny, nz);
        end
    end
end
