function [tau, tau_f] = ID_flex(model, q, qd, qdd, eta, eta_d, eta_dd, f_ext, f_modal)

NB = model.NB;
fixed = fixed_joint_mask(model);

if nargin < 8 || isempty(f_ext)
    f_ext = cell(NB,1); 
end
if nargin < 9 || isempty(f_modal)
    f_modal = cell(NB,1); 
end

f_ext(end+1:NB) = {[]}; 
f_modal(end+1:NB) = {[]};

gravity = zeros(6,1);

if isfield(model,'gravity')
    gravity = model.gravity(:); 
end
if numel(gravity)==3
    gravity = [zeros(3,1);gravity]; 
end

for i=1:NB
    p=model.parent(i);

    if is_flexible_link(model,i)
     % TODO split B{i} and eta_i,eta_d_i into 2 functions
    % Change function name 'local_modal_state'
        %[bodyPhi{i},e{i},ed{i}] = local_modal_state(model,i,eta,eta_d);
        bodyPhi{i} = model.Phi{i};
        n = size(bodyPhi{i}, 2);
        e{i}=local_link_vector(model,i,eta,n);
        ed{i}=local_link_vector(model,i,eta_d,n);
        edd{i}=local_link_vector(model,i,eta_dd,n);
    else
        bodyPhi{i}=zeros(6,0); 
        e{i}=zeros(0,1); 
        ed{i}=e{i}; 
        edd{i}=e{i};
    end
    if fixed(i)
        [XJ,~]=jcalc(model.jtype{i},0); 
        S{i}=zeros(6,1);
    else
        [XJ,S{i}]=jcalc(model.jtype{i},q(i));
    end

    Phi{i}=zeros(6,0); 
    Xdef=eye(6);
    if p~=0 && is_flexible_link(model,p)
        % Change Name
        Phi{i}=endpoint_shape(model,i,p,bodyPhi{p},numel(e{p}));
        Xdef=fcalc(Phi{i},e{p});
    end
    Xup{i}=XJ*Xdef*model.Xtree{i};
    A{i}=XJ*spatial_rotation(Xdef*model.Xtree{i});
    vF=zeros(6,1); 
    aF=zeros(6,1); 
    vFpre=zeros(6,1);
    if p==0
        vp=zeros(6,1); 
        ap=-gravity; 
        X0{i}=Xup{i};
    else
        vp=v{p}; 
        ap=a{p}; 
        X0{i}=Xup{i}*X0{p};
        if is_flexible_link(model,p)
            vFpre=spatial_rotation(Xdef*model.Xtree{i})*Phi{i}*ed{p};
            vF=A{i}*Phi{i}*ed{p};
            aF=A{i}*Phi{i}*edd{p};
        end
    end
    vRel=S{i}*qd(i)+vF;
    v{i}=Xup{i}*vp+vRel;
    % 
    gamma=flex_connection_bias(v{i},vRel,vF,XJ,vFpre);
    a{i}=Xup{i}*ap+S{i}*qdd(i)+aF+gamma;
    if is_flexible_link(model,i)
        % Cancel function
        body=local_body(model,i,true);
        % Classes function
        [Mrr,Mrf,Mff,cr,cf]=body.update_matrices(v{i},e{i},ed{i});

        f{i}=Mrr*a{i}+Mrf*edd{i}+cr;
        tau_f{i}=Mrf.'*a{i}+Mff*edd{i}+cf;
    else
        body=local_body(model,i,false); 
        Mrr=body.get_spatial_inertia();
        f{i}=Mrr*a{i}+crf(v{i})*Mrr*v{i}; 
        tau_f{i}=zeros(0,1);
    end
    if ~isempty(f_ext{i})
        % External wrenches are body-frame by default.
        f{i}=f{i}-f_ext{i}(:);
        tau_f{i}=tau_f{i}-bodyPhi{i}.'*f_ext{i}(:);
    end
    if ~isempty(f_modal{i})
        tau_f{i}=tau_f{i}-f_modal{i}(:);
    end
end
tau=0*qdd; 
for i=NB:-1:1
    tau(i)=S{i}.'*f{i};
    p=model.parent(i);
    if p~=0
        f{p}=f{p}+Xup{i}.'*f{i};
        if is_flexible_link(model,p)
            tau_f{p}=tau_f{p}+Phi{i}.'*A{i}.'*f{i};
        end
    end
end
end

