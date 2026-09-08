function R6 = spatial_rotation(X)
R = X(1:3, 1:3);
R6 = blkdiag(R, R);
end