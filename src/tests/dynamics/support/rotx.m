function X = rotx(q)
% Passive spatial rotation; q is measured in radians.
E=get_rotation_matrix([1;0;0],-q);
X=blkdiag(E,E);
end
