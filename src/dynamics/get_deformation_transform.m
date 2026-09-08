function [R, T, delta] = get_deformation_transform(Phi, eta)
%GET_DEFORMATION_TRANSFORM
% 将模态小变形 delta = Phi*eta 转换为旋转矩阵和齐次变换矩阵。
%
% 输入：
%   Phi   : 6×n 模态矩阵或变形映射矩阵
%   eta   : n×1 模态坐标

%
% 输出：
%   R     : 3×3 旋转矩阵
%   T     : 4×4 齐次变换矩阵
%   delta : 6×1 小变形向量 Phi*eta


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
    angle = norm(theta);

    if angle < 1.0e-8
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

    T = [R,       u;
         0, 0, 0, 1];

end
