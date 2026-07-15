function [R, T, delta] = get_deformation_transform(Phi, eta, method)
%GET_DEFORMATION_TRANSFORM
% 将模态小变形 delta = Phi*eta 转换为旋转矩阵和齐次变换矩阵。
%
% 输入：
%   Phi   : 6×n 模态矩阵或变形映射矩阵
%   eta   : n×1 模态坐标
%   method:
%       'exp'    使用旋转向量指数映射，默认选项
%       'linear' 使用小转角一阶近似
%
% 输出：
%   R     : 3×3 旋转矩阵
%   T     : 4×4 齐次变换矩阵
%   delta : 6×1 小变形向量 Phi*eta
%
% 小变形排列约定：
%
%       delta = [theta_x;
%                theta_y;
%                theta_z;
%                u_x;
%                u_y;
%                u_z]
%
% 齐次变换的作用约定：
%
%       p_deformed = R*p_reference + u

    if nargin < 3 || isempty(method)
        method = 'exp';
    end


    if ~isnumeric(Phi) || ~ismatrix(Phi)
        error('Phi 必须是数值矩阵。');
    end

    if ~isnumeric(eta) || ~isvector(eta)
        error('eta 必须是数值向量。');
    end

    eta = eta(:);

    if size(Phi, 2) ~= length(eta)
        error(['Phi 的列数必须等于 eta 的长度。', ...
               '当前 size(Phi,2)=%d，length(eta)=%d。'], ...
               size(Phi, 2), length(eta));
    end

    delta = Phi * eta;

    if length(delta) ~= 6
        error(['Phi*eta 必须是 6×1 向量，排列为 ', ...
               '[theta_x;theta_y;theta_z;u_x;u_y;u_z]。']);
    end

    theta = delta(1:3);
    u     = delta(4:6);


    theta_tilde = get_skew_sym_matrix(theta);


    switch lower(method)

        case 'exp'
            % R = exp([theta]x)，使用 Rodrigues 公式

            angle = norm(theta);

            if angle < 1.0e-8
                % 泰勒展开，避免 angle 接近零时出现数值问题
                A = 1 ...
                    - angle^2 / 6 ...
                    + angle^4 / 120;

                B = 1/2 ...
                    - angle^2 / 24 ...
                    + angle^4 / 720;
            else
                A = sin(angle) / angle;
                B = (1 - cos(angle)) / angle^2;
            end

            R = eye(3) ...
                + A * theta_tilde ...
                + B * theta_tilde^2;

        case 'linear'
            % 小转角一阶近似
            %
            % 注意：该矩阵一般不严格满足 R'*R = I
            R = eye(3) + theta_tilde;

        otherwise
            error('method 只能取 ''exp'' 或 ''linear''。');
    end


    T = [R,       u;
         0, 0, 0, 1];

end


function S = get_skew_sym_matrix(v)
%GET_SKEW_SYM_MATRIX 返回三维向量对应的反对称矩阵

    if ~isnumeric(v) || numel(v) ~= 3
        error('输入 v 必须包含三个数值元素。');
    end

    v = v(:);

    S = [ 0,    -v(3),  v(2);
          v(3),  0,    -v(1);
         -v(2),  v(1),  0    ];

end