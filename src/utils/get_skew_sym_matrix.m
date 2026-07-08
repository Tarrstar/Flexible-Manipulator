function S = get_skew_sym_matrix(v)
    if ~isvector(v) || numel(v) ~= 3
        error('the input has to be a 3d vector');
    end

    v = v(:);
    
    S = [  0,    -v(3),  v(2);
          v(3),  0,    -v(1);
         -v(2),  v(1),  0   ];
end