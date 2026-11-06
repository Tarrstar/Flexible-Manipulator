function [dtau_dq,dtau_dqd] = ID_flex_FO(model,q,qd,qdd,eta,eta_d,eta_dd,f_ext,f_modal)
%ID_FLEX_FO Composite first derivatives with joint/modal blocks assembled together.
% Rows [tau;tau_f], columns [q;eta] or [qd;eta_d]. The eight mathematical
% blocks share matrix products and writes; they are submatrices of these outputs.
% modal{i} includes the NB offset, cols{i} selects local connection sources,
% rows{i} selects subtree force rows, and ancestor{i} selects strict ancestors.

NB=model.NB; q=q(:); qd=qd(:); qdd=qdd(:);
fixed=fixed_joint_mask(model);
if nargin<8 || isempty(f_ext), f_ext=cell(NB,1); end
f_ext(end+1:NB)={[]};
inertial=isfield(model,'f_ext_frame') && strcmpi(model.f_ext_frame,'inertial');
gravity=zeros(6,1);
if isfield(model,'gravity'), gravity=model.gravity(:); end
if numel(gravity)==3, gravity=[zeros(3,1);gravity]; end

flexible=false(NB,1); 
nf=0;
for i=1:NB
    flexible(i)=is_flexible_link(model,i);
    body{i}=local_body(model,i,flexible(i));
    if flexible(i)
        bodyPhi{i}=local_modal_state(model,i); 
        n=body{i}.nEta;
        e{i}=local_link_vector(model,i,eta,n);
        ed{i}=local_link_vector(model,i,eta_d,n);
        edd{i}=local_link_vector(model,i,eta_dd,n);
    else
        bodyPhi{i}=zeros(6,0); 
        e{i}=zeros(0,1); 
        ed{i}=e{i}; 
        edd{i}=e{i};
    end
    modal{i}=NB+nf+(1:numel(e{i}));
    nf=nf+numel(e{i});
end
nv=NB+nf;
J=cell(NB,1); ancestor=J; rows=J; cols=J;

for i=1:NB
    p=model.parent(i); parentFlexible=p~=0 && flexible(p);
    if fixed(i)
        [XJ,~]=jcalc(model.jtype{i},0); 
        Si=zeros(6,1);
    else
        [XJ,Si]=jcalc(model.jtype{i},q(i));
    end
    S_cross=crm(Si); 
    modes=[]; 
    Phi{i}=zeros(6,0);
    if parentFlexible
        modes=modal{p};
        Phi{i}=endpoint_shape(model,i,p,bodyPhi{p},numel(e{p}));
        [Xdef,Xdef_qe,gdef_qe]=deformation_terms(Phi{i},e{p});
    else
        Xdef=eye(6);
    end
    cols{i}=[i,modes];
    m=numel(cols{i});
    Xup=XJ*Xdef*model.Xtree{i};
    Rpre=spatial_rotation(Xdef*model.Xtree{i}); 
    A{i}=XJ*Rpre;
    if p==0
        Xup0{i}=Xup; 
        vp=zeros(6,1); 
        ap=-gravity;
        J{i}=zeros(6,nv);
        ancestor{i}=[];
    else
        Xup0{i}=Xup*Xup0{p}; 
        vp=v{p}; 
        ap=a{p};
        J{i}=J{p};
        ancestor{i}=[ancestor{p},p,modal{p}];
    end
    X=Xup0{i}; 
    T=inv(X); 
    S{i}=T*Si;
    g{i}=zeros(6,m); 
    g{i}(:,1)=S{i};
    if parentFlexible
        g{i}(:,2:end)=-T*XJ*gdef_qe; 
    end
    F=A{i}*Phi{i}; 
    Fpre=Rpre*Phi{i}; 
    E{i}=T*F;
    V{i}=[S{i},E{i}]; 
    J{i}(:,cols{i})=J{i}(:,cols{i})+V{i};
    rows{i}=[i,modal{i}];
    if parentFlexible
        vF=F*ed{p}; 
        aF=F*edd{p}; 
        vpre=Fpre*ed{p};
    else
        vF=zeros(6,1);
        aF=vF; 
        vpre=vF;
    end
    vF_q=zeros(6,m); 
    aF_q=vF_q; 
    vpre_q=vF_q;
    vF_q(:,1)=-S_cross*vF; 
    aF_q(:,1)=-S_cross*aF;
    A_q{i}=-S_cross*A{i}; 
    A_qe{i}=zeros(6,6,numel(modes));
    for j=1:numel(modes)
        dR=spatial_rotation(Xdef_qe(:,:,j)*model.Xtree{i});
        dA=XJ*dR; 
        A_qe{i}(:,:,j)=dA;
        vF_q(:,j+1)=dA*Phi{i}*ed{p};
        aF_q(:,j+1)=dA*Phi{i}*edd{p};
        vpre_q(:,j+1)=dR*Phi{i}*ed{p};
    end

    % Only local edge-source derivatives, not all ancestor columns.
    vJ_body=Si*qd(i); 
    v_J=S{i}*qd(i)+T*vF; 
    a_J=S{i}*qdd(i)+T*aF;
    v_J_q=-crm(v_J)*g{i}+T*vF_q;
    a_J_q=-crm(a_J)*g{i}+T*aF_q;
    W=get_skew_sym_matrix(vpre(1:3));
    preBias=[zeros(3,1);W*vpre(4:6)];
    preBias_v=[zeros(3,6);-get_skew_sym_matrix(vpre(4:6)),W];
    kappa_b=-crm(vJ_body)*vF-XJ*preBias; 
    kappa=T*kappa_b;
    kappa_q=-crm(vJ_body)*vF_q-XJ*preBias_v*vpre_q;
    kappa_q(:,1)=kappa_q(:,1)+S_cross*XJ*preBias;
    kappa_q=-crm(kappa)*g{i}+T*kappa_q;
    kappa_qd=zeros(6,m); 
    kappa_qd(:,1)=-S_cross*vF;
    if parentFlexible
        kappa_qd(:,2:end)=-crm(vJ_body)*F-XJ*preBias_v*Fpre;
    end
    v{i}=vp+v_J; 
    a{i}=ap+a_J+crm(v{i})*v_J+kappa;
    P{i}=v_J_q+crm(v{i})*g{i};
    H{i}=a_J_q+crm(vp)*v_J_q+kappa_q+crm(a{i})*g{i}+crm(v{i})*P{i};
    Hv{i}=crm(vp+v{i})*V{i}+T*kappa_qd;

    vb=X*v{i}; 
    ab=X*a{i};
    if flexible(i)
        [Mrr,Mrf,C,C_v,C_qe,C_qed,c_v,c_qe,c_qed,Mrr_qe,Mrf_qe]= ...
            body_terms(body{i},vb,e{i},ed{i});
        fb=Mrr*ab+Mrf*edd{i}+C;
        intrinsic_q{i}=zeros(6,numel(e{i}));
        intrinsic_f_q{i}=c_qe;
        for j=1:numel(e{i})
            intrinsic_q{i}(:,j)=Mrr_qe(:,:,j)*ab+Mrf_qe(:,:,j)*edd{i}+C_qe(:,j);
            intrinsic_f_q{i}(:,j)=intrinsic_f_q{i}(:,j)+Mrf_qe(:,:,j).'*ab;
        end
        intrinsic_q{i}=X.'*intrinsic_q{i}; 
        intrinsic_qd{i}=X.'*C_qed;
        intrinsic_f_qd{i}=c_qed;
        modalIC{i}=Mrf.'*X;
        modalBC{i}=c_v*X-modalIC{i}*crm(v{i});
        modalGC{i}=zeros(numel(e{i}),6);
    else
        Mrr=body{i}.get_spatial_inertia(); 
        h=Mrr*vb;
        C_v=crf(vb)*Mrr+force_right(h); 
        fb=Mrr*ab+crf(vb)*h;
        modalIC{i}=zeros(0,6); 
        modalBC{i}=modalIC{i}; 
        modalGC{i}=modalIC{i};
    end
    IC{i}=X.'*Mrr*X;
    BC{i}=X.'*C_v*X-IC{i}*crm(v{i});
    f_transport{i}=X.'*fb;
    if ~isempty(f_ext{i})
        wrench=f_ext{i}(:);
        if inertial
            wrench=X.'\wrench;
            if flexible(i)
                modalGC{i}=bodyPhi{i}.'*force_right(wrench)*X;
            end
        end
        fb=fb-wrench;
    end
    f{i}=X.'*fb;
    if ~inertial
        f_transport{i}=f{i};
    end
end

% Composite inertias first; assemble low-rank derivative blocks afterward.
dtau_dq=zeros(nv); 
dtau_dqd=dtau_dq;
tmp1=zeros(nv,6); 
tmp2=tmp1; 
tmp3=tmp1;
for i=NB:-1:1
    p=model.parent(i); 
    rr=[i,modal{i}];
    transport=force_right(f_transport{i}-f{i});
    tmp1(rr,:)=[S{i}.'*IC{i};modalIC{i}];
    tmp2(rr,:)=[S{i}.'*BC{i};modalBC{i}];
    tmp3(rr,:)=[S{i}.'*transport;modalGC{i}];
    r=rows{i}; 
    c=cols{i};
    dtau_dq(r,c)=dtau_dq(r,c)+tmp2(r,:)*P{i}+tmp1(r,:)*H{i}+tmp3(r,:)*g{i};
    dtau_dqd(r,c)=dtau_dqd(r,c)+tmp2(r,:)*V{i}+tmp1(r,:)*Hv{i};
    r=ancestor{i};
    q_tmp=BC{i}*P{i}+IC{i}*H{i}+force_right(f_transport{i})*g{i};
    qd_tmp=BC{i}*V{i}+IC{i}*Hv{i};
    dtau_dq(r,c)=dtau_dq(r,c)+J{i}(:,r).'*q_tmp;
    dtau_dqd(r,c)=dtau_dqd(r,c)+J{i}(:,r).'*qd_tmp;
    if p~=0 && flexible(p)
        % The parent's modal projection at this edge also changes locally.
        dtau_dq(modal{p},c)=dtau_dq(modal{p},c)-E{i}.'*force_right(f{i})*g{i};
        fb=Xup0{i}.'\f{i};
        dtau_dq(modal{p},i)=dtau_dq(modal{p},i)+Phi{i}.'*A_q{i}.'*fb;
        for j=1:numel(modal{p})
            col=modal{p}(j);
            dtau_dq(modal{p},col)=dtau_dq(modal{p},col)+Phi{i}.'*A_qe{i}(:,:,j).'*fb;
        end
    end
    if flexible(i)
        % Intrinsic inertia/bias sources act only at this body, not on children.
        r=[i,ancestor{i}]; 
        c=modal{i};
        dtau_dq(r,c)=dtau_dq(r,c)+J{i}(:,r).'*intrinsic_q{i};
        dtau_dqd(r,c)=dtau_dqd(r,c)+J{i}(:,r).'*intrinsic_qd{i};
        dtau_dq(c,c)=dtau_dq(c,c)+intrinsic_f_q{i};
        dtau_dqd(c,c)=dtau_dqd(c,c)+intrinsic_f_qd{i};
    end
    if p~=0
        if flexible(p)
            modalIC{p}=modalIC{p}+E{i}.'*IC{i};
            modalBC{p}=modalBC{p}+E{i}.'*BC{i};
            modalGC{p}=modalGC{p}+E{i}.'*transport;
        end
        IC{p}=IC{p}+IC{i}; 
        BC{p}=BC{p}+BC{i};
        f{p}=f{p}+f{i}; 
        f_transport{p}=f_transport{p}+f_transport{i};
        rows{p}=[rows{p},rows{i}];
    end
end
end
