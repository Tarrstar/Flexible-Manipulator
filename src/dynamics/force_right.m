function B=force_right(h)
% B*dv = crf(dv)*h, avoiding six basis-vector cross-matrix constructions.
ha=get_skew_sym_matrix(h(1:3)); hl=get_skew_sym_matrix(h(4:6));
B=[-ha,-hl;-hl,zeros(3)];
end

