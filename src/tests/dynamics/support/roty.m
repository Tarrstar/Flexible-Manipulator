function X = roty(q)
% Passive spatial rotation; q is measured in radians.
E=get_rotation_matrix([0;1;0],-q);
X=blkdiag(E,E);
end
