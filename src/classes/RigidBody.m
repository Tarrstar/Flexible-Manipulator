classdef RigidBody < handle
    %RIGIDBODY Rigid body with directly specified principal inertias.
    %
    % params.mass       scalar mass
    % params.inertia    3-vector [Ixx Iyy Izz] or a 3-by-3 inertia matrix
    % params.com        optional COM position in the body frame (default 0)

    properties (SetAccess = private)
        mass
        inertia
        com
    end

    methods
        function obj = RigidBody(params)
            if nargin == 0
                return
            end
            if ~isstruct(params) || ~isfield(params, 'mass')
                error('RigidBody requires params.mass.');
            end
            if ~isscalar(params.mass) || ~isreal(params.mass) || ...
                    ~isfinite(params.mass) || params.mass <= 0
                error('params.mass must be a positive finite scalar.');
            end
            if isfield(params, 'inertia')
                I = params.inertia;
            elseif all(isfield(params, {'Ixx', 'Iyy', 'Izz'}))
                I = [params.Ixx; params.Iyy; params.Izz];
            else
                error('Provide params.inertia or params.Ixx/Iyy/Izz.');
            end
            if isvector(I) && numel(I) == 3
                I = diag(I(:));
            end
            if ~isequal(size(I), [3, 3]) || ~isreal(I) || any(~isfinite(I(:)))
                error('params.inertia must be a finite 3-vector or 3-by-3 matrix.');
            end
            obj.mass = params.mass;
            obj.inertia = 0.5 * (I + I.');
            if isfield(params, 'com') && ~isempty(params.com)
                obj.com = params.com(:);
            else
                obj.com = zeros(3, 1);
            end
            if numel(obj.com) ~= 3 || any(~isfinite(obj.com))
                error('params.com must be a finite 3-vector.');
            end
        end

        function I = get_spatial_inertia(obj)
            % Spatial ordering is [angular velocity; linear velocity].
            c = get_skew_sym_matrix(obj.com);
            Iorigin = obj.inertia - obj.mass * c * c;
            I = [Iorigin, obj.mass * c;
                 -obj.mass * c, obj.mass * eye(3)];
        end
    end
end
