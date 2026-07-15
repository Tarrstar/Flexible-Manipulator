function  tau = ID_flex( model, q, qd, qdd, eta, eat_d, f_ext )

a_grav = get_gravity(model);

for i = 1:model.NB
  [ XJ, S{i}] = jcalc( model.jtype{i}, q(i) );
  Phi{i} = model.Phi{i};
  [R{i}, T{i}, delta{i}]  = get_deformation_transformarion(Phi{i}*eta);
  Xup{i} = T{i}* XJ * model.Xtree{i};
 
  vJ = S{i}*qd(i);
  vF = R{i}*Phi{i-1}*eta_d(i-1);
  
  if model.parent(i) == 0
    v{i} = vJ+vF;
    a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i)+R{i}*Phi{i-1}*eta_dd(i-1);
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ+vF;
    a{i} = Xup{i}*a{model.parent(i)} + S{i}*qdd(i) + crm(v{i})*(vJ+vF);
  end
  model.I_rr{i} = ;
  model.I_re{i} = ;
  model.M_ee{i} = ;
  mode.c{i} = ;
  f{i} = model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i};
end

if nargin == 5
  f = apply_external_forces( model.parent, Xup, f, f_ext );
end

for i = model.NB:-1:1
  if model.link_type == 'flexible'
        K_plus{i} = I{i};
        beta{i} = C;
  elseif model.link_type == 'rigid'
        K_plus{i} = I{i};
        beta{i} = C{i};
  end
  
  tau(i,1) = S{i}' * f{i};
  if model.parent(i) ~= 0
    f{model.parent(i)} = f{model.parent(i)} + Xup{i}'*f{i};
  end
end
