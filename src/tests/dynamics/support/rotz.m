function X = rotz(q)
% Spatial z-axis rotation transform used by jcalc in the dynamics test.
E = get_rotation_matrix([0; 0; 1], -q);
X = blkdiag(E, E);
end
