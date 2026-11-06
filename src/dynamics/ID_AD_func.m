function ID_d = ID_AD_func(model, f_ext)

if nargin<2
    f_ext=[]; 
end

x=casadi.MX.sym('position',model.ndof);
v=casadi.MX.sym('velocity',model.ndof);
a=casadi.MX.sym('acceleration',model.ndof);
q=cell(model.NB,1); qd=q; qdd=q;
for i=1:model.NB
    k=find(model.active_joints==i,1);
    if isempty(k)
        q{i}=0*x(1); qd{i}=0*v(1); qdd{i}=0*a(1);
    else
        q{i}=x(k); qd{i}=v(k); qdd{i}=a(k);
    end
end
m=model.joint_num + (1:model.modal_num);
[tau,tf]=ID_flex(model,vertcat(q{:}),vertcat(qd{:}),vertcat(qdd{:}), ...
    x(m),v(m),a(m),f_ext);
value=[tau(model.active_joints);vertcat(tf{:})];
ID_d=casadi.Function('flex_inverse_dynamics',{x,v,a}, ...
    {value,jacobian(value,x),jacobian(value,v),jacobian(value,a)}, ...
    {'q','v','a'},{'force','Jq','Jv','Ja'});
end
