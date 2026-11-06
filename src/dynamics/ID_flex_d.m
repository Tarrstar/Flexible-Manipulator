function [tau_q, tau_qd] = ID_flex_d(model,q,qd,qdd,eta,eta_d,eta_dd,f_ext,f_modal)
%ID_FLEX_D Analytic inverse-dynamics position and velocity Jacobians.
% Rows [joint torques; modes]; columns [q; eta] and [qd; eta_d].
% Inputs are assumed valid. f_ext is a cell array of constant wrenches.
% Uses the default floating-frame kernel. Constant f_modal is accepted for
% compatibility with cID_flex but has no contribution to these derivatives.

NB=model.NB; q=q(:); qd=qd(:); qdd=qdd(:);
fixed=fixed_joint_mask(model);
if nargin<8 || isempty(f_ext), f_ext=cell(NB,1); end
f_ext(end+1:NB)={[]};
inertial=isfield(model,'f_ext_frame') && strcmpi(model.f_ext_frame,'inertial');
gravity=zeros(6,1);
if isfield(model,'gravity'), gravity=model.gravity(:); end
if numel(gravity)==3, gravity=[zeros(3,1);gravity]; end

flexible=false(NB,1); nf=0;

for i=1:NB
    flexible(i)=is_flexible_link(model,i);
    body{i}=local_body(model,i,flexible(i));
    if flexible(i)
        bodyPhi{i}=local_modal_state(model,i);
        e{i}=local_link_vector(model,i,eta,body{i}.nEta);
        ed{i}=local_link_vector(model,i,eta_d,body{i}.nEta);
        edd{i}=local_link_vector(model,i,eta_dd,numel(e{i}));
    else
        bodyPhi{i}=zeros(6,0); e{i}=zeros(0,1); ed{i}=e{i}; edd{i}=e{i};
    end
    modal{i}=nf+(1:numel(e{i})); nf=nf+numel(e{i});
end

for i=1:NB
    p=model.parent(i); parentFlexible=p~=0 && flexible(p);
    if fixed(i)
        [XJ,~]=jcalc(model.jtype{i},0); S{i}=zeros(6,1);
    else
        [XJ,S{i}]=jcalc(model.jtype{i},q(i));
    end
    S_cross=crm(S{i});
    XJ_q=-S_cross*XJ;
    modes=[];
    if parentFlexible
        modes=modal{p};
        Phi{i}=endpoint_shape(model,i,p,bodyPhi{p},numel(e{p}));
        [Xdef,Xdef_qe,gdef_qe]=deformation_terms(Phi{i},e{p});
    else
        Xdef=eye(6);
    end
    Xup{i}=XJ*Xdef*model.Xtree{i};
    Xup_q{i}=-S_cross*Xup{i};
    Rpre=spatial_rotation(Xdef*model.Xtree{i});
    A{i}=XJ*Rpre; A_q{i}=-S_cross*A{i};
    Xup_qe{i}=zeros(6,6,numel(modes)); A_qe{i}=Xup_qe{i};

    if p~=0
        vp=v{p}; ap=a{p};
        v_q{i}=Xup{i}*v_q{p}; v_qd{i}=Xup{i}*v_qd{p};
        v_qe{i}=Xup{i}*v_qe{p}; v_qed{i}=Xup{i}*v_qed{p};
        a_q{i}=Xup{i}*a_q{p}; a_qd{i}=Xup{i}*a_qd{p};
        a_qe{i}=Xup{i}*a_qe{p}; a_qed{i}=Xup{i}*a_qed{p};
    else
        vp=zeros(6,1); ap=-gravity;
        v_q{i}=zeros(6,NB); v_qd{i}=v_q{i};
        v_qe{i}=zeros(6,nf); v_qed{i}=v_qe{i};
        a_q{i}=zeros(6,NB); a_qd{i}=a_q{i};
        a_qe{i}=zeros(6,nf); a_qed{i}=a_qe{i};
    end
    if parentFlexible
        Fpre=Rpre*Phi{i}; F=A{i}*Phi{i};
        vpre=Fpre*ed{p}; vF=F*ed{p}; aF=F*edd{p};
    else
        vpre=zeros(6,1); vF=vpre; aF=vpre;
    end
    vJ=S{i}*qd(i); v_rel=vJ+vF;
    v{i}=Xup{i}*vp+v_rel;
    V_cross=crm(v{i}); J_cross=crm(vJ);
    Wpre=get_skew_sym_matrix(vpre(1:3));
    % d(w x u) = -[u]x dw + [w]x du; w x w is identically zero.
    preBias=[zeros(3,1);Wpre*vpre(4:6)];
    preBias_v=[zeros(3,6);-get_skew_sym_matrix(vpre(4:6)),Wpre];
    a{i}=Xup{i}*ap+S{i}*qdd(i)+aF+V_cross*v_rel-J_cross*vF-XJ*preBias;

    vF_q=-S_cross*vF;
    v_q{i}(:,i)=v_q{i}(:,i)+Xup_q{i}*vp+vF_q;
    v_qd{i}(:,i)=v_qd{i}(:,i)+S{i};
    if parentFlexible, v_qed{i}(:,modes)=v_qed{i}(:,modes)+F; end
    % Cross antisymmetry turns all-column crm(dv)*rel into -crm(rel)*dv.
    v_rel_cross=-crm(v_rel);
    a_q{i}=a_q{i}+v_rel_cross*v_q{i};
    a_qd{i}=a_qd{i}+v_rel_cross*v_qd{i};
    a_q{i}(:,i)=a_q{i}(:,i)+Xup_q{i}*ap-S_cross*aF ...
        +(V_cross-J_cross)*vF_q-XJ_q*preBias;
    a_qd{i}(:,i)=a_qd{i}(:,i)+V_cross*S{i}-S_cross*vF;

    for j=1:numel(modes)
        col=modes(j);
        dX=XJ*Xdef_qe(:,:,j)*model.Xtree{i};
        dR=spatial_rotation(Xdef_qe(:,:,j)*model.Xtree{i});
        dA=XJ*dR;
        Xup_qe{i}(:,:,j)=dX; A_qe{i}(:,:,j)=dA;
        dvpre=dR*Phi{i}*ed{p}; dvF=dA*Phi{i}*ed{p};
        v_qe{i}(:,col)=v_qe{i}(:,col)+dX*vp+dvF;
        a_qe{i}(:,col)=a_qe{i}(:,col)+dX*ap+dA*Phi{i}*edd{p} ...
            +(V_cross-J_cross)*dvF-XJ*preBias_v*dvpre;
    end
    a_qe{i}=a_qe{i}+v_rel_cross*v_qe{i};
    a_qed{i}=a_qed{i}+v_rel_cross*v_qed{i};
    if parentFlexible
        a_qed{i}(:,modes)=a_qed{i}(:,modes) ...
            +(V_cross-J_cross)*F-XJ*preBias_v*Fpre;
    end

    if flexible(i)
        [Mrr,Mrf,C,C_v,C_qe,C_qed,c_v,c_qe,c_qed,Mrr_qe,Mrf_qe]= ...
            body_terms(body{i},v{i},e{i},ed{i});
        f{i}=Mrr*a{i}+Mrf*edd{i}+C;
        f_q{i}=Mrr*a_q{i}+C_v*v_q{i};
        f_qd{i}=Mrr*a_qd{i}+C_v*v_qd{i};
        f_qe{i}=Mrr*a_qe{i}+C_v*v_qe{i};
        f_qed{i}=Mrr*a_qed{i}+C_v*v_qed{i};
        MrfT=Mrf.';
        tau_f_q{i}=MrfT*a_q{i}+c_v*v_q{i};
        tau_f_qd{i}=MrfT*a_qd{i}+c_v*v_qd{i};
        tau_f_qe{i}=MrfT*a_qe{i}+c_v*v_qe{i};
        tau_f_qed{i}=MrfT*a_qed{i}+c_v*v_qed{i};
        for j=1:numel(e{i})
            col=modal{i}(j);
            f_qe{i}(:,col)=f_qe{i}(:,col)+Mrr_qe(:,:,j)*a{i} ...
                +Mrf_qe(:,:,j)*edd{i}+C_qe(:,j);
            tau_f_qe{i}(:,col)=tau_f_qe{i}(:,col)+Mrf_qe(:,:,j).'*a{i}+c_qe(:,j);
        end
        f_qed{i}(:,modal{i})=f_qed{i}(:,modal{i})+C_qed;
        tau_f_qed{i}(:,modal{i})=tau_f_qed{i}(:,modal{i})+c_qed;
    else
        Mrr=body{i}.get_spatial_inertia();
        h=Mrr*v{i}; V_star=crf(v{i});
        C_v=V_star*Mrr+force_right(h);
        f{i}=Mrr*a{i}+V_star*h;
        f_q{i}=Mrr*a_q{i}+C_v*v_q{i};
        f_qd{i}=Mrr*a_qd{i}+C_v*v_qd{i};
        f_qe{i}=Mrr*a_qe{i}+C_v*v_qe{i};
        f_qed{i}=Mrr*a_qed{i}+C_v*v_qed{i};
        tau_f_q{i}=zeros(0,NB); tau_f_qd{i}=tau_f_q{i};
        tau_f_qe{i}=zeros(0,nf); tau_f_qed{i}=tau_f_qe{i};
    end

    if inertial
        % g0 represents dX0=crm(g0)*X0. Six rows suffice instead of
        % 6-by-6-by-global-coordinate tensors; dw=icrf(w)*g0.
        g0_q{i}=zeros(6,NB); g0_qe{i}=zeros(6,nf);
        X0{i}=Xup{i};
        if p~=0
            X0{i}=Xup{i}*X0{p};
            g0_q{i}=Xup{i}*g0_q{p}; g0_qe{i}=Xup{i}*g0_qe{p};
        end
        g0_q{i}(:,i)=g0_q{i}(:,i)-S{i};
        if parentFlexible, g0_qe{i}(:,modes)=g0_qe{i}(:,modes)+XJ*gdef_qe; end
    end
    if ~isempty(f_ext{i})
        wrench=f_ext{i}(:);
        if inertial
            wrench=X0{i}.'\wrench;
            wrench_cross=force_right(wrench);
            wrench_q=wrench_cross*g0_q{i}; wrench_qe=wrench_cross*g0_qe{i};
            f_q{i}=f_q{i}-wrench_q; f_qe{i}=f_qe{i}-wrench_qe;
            tau_f_q{i}=tau_f_q{i}-bodyPhi{i}.'*wrench_q;
            tau_f_qe{i}=tau_f_qe{i}-bodyPhi{i}.'*wrench_qe;
        end
        f{i}=f{i}-wrench;
    end
end

tau_q=zeros(NB,NB); tau_qd=tau_q; tau_qe=zeros(NB,nf); tau_qed=tau_qe;
for i=NB:-1:1
    St=S{i}.';
    tau_q(i,:)=St*f_q{i}; tau_qd(i,:)=St*f_qd{i};
    tau_qe(i,:)=St*f_qe{i}; tau_qed(i,:)=St*f_qed{i};
    p=model.parent(i);
    if p~=0
        XT=Xup{i}.';
        f_q{p}=f_q{p}+XT*f_q{i}; f_qd{p}=f_qd{p}+XT*f_qd{i};
        f_qe{p}=f_qe{p}+XT*f_qe{i}; f_qed{p}=f_qed{p}+XT*f_qed{i};
        f_q{p}(:,i)=f_q{p}(:,i)+Xup_q{i}.'*f{i};
        if flexible(p)
            PhiT=Phi{i}.'; projection=PhiT*A{i}.';
            tau_f_q{p}=tau_f_q{p}+projection*f_q{i};
            tau_f_qd{p}=tau_f_qd{p}+projection*f_qd{i};
            tau_f_qe{p}=tau_f_qe{p}+projection*f_qe{i};
            tau_f_qed{p}=tau_f_qed{p}+projection*f_qed{i};
            tau_f_q{p}(:,i)=tau_f_q{p}(:,i)+PhiT*A_q{i}.'*f{i};
            for j=1:numel(modal{p})
                col=modal{p}(j);
                f_qe{p}(:,col)=f_qe{p}(:,col)+Xup_qe{i}(:,:,j).'*f{i};
                tau_f_qe{p}(:,col)=tau_f_qe{p}(:,col)+PhiT*A_qe{i}(:,:,j).'*f{i};
            end
        end
        f{p}=f{p}+XT*f{i};
    end
end
tau_q=[tau_q,tau_qe;vertcat(tau_f_q{:}),vertcat(tau_f_qe{:})];
tau_qd=[tau_qd,tau_qed;vertcat(tau_f_qd{:}),vertcat(tau_f_qed{:})];
end

