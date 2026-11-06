function [M,H,C,Cv,Ceta,Ced,cv,ceta,ced,dM,dH]=body_terms(body,v,eta,ed)
% Shared cached contractions for primal inertia and both derivative classes.
n=body.nEta; qa=[1;eta]; qad=[0;ed]; T=zeros(3); Td=T;
Te=zeros(3,3,n); Tde=Te; P=cell(3); V=P;
W=get_skew_sym_matrix(v(1:3)); W2=W*W;
cv=zeros(n,6); cv(:,4:6)=body.S_bar.'*W;
ceta=body.K_ff; ced=body.D_ff;
for k=1:3
    for l=1:3
        product=body.shapeProducts{k,l}; pq=product*qa; pqd=product*qad;
        T(k,l)=qa.'*pq; Td(k,l)=qad.'*pq+qa.'*pqd;
        Te(k,l,:)=reshape(pq(2:end)+product(:,2:end).'*qa,1,1,n);
        Tde(k,l,:)=reshape(pqd(2:end)+product(:,2:end).'*qad,1,1,n);
        P{k,l}=reshape(pq(2:end),n,1); V{k,l}=reshape(pqd(2:end),n,1);
        ceta=ceta+W2(k,l)*product(2:end,2:end);
        ced=ced+2*W(k,l)*product(2:end,2:end);
    end
end
J=inertia_from_moment(T); Jd=inertia_from_moment(Td);
Hm=[(P{3,2}-P{2,3}).';(P{1,3}-P{3,1}).';(P{2,1}-P{1,2}).'];
st=get_skew_sym_matrix(body.S_t0+body.S_bar*eta);
sdt=get_skew_sym_matrix(body.S_bar*ed);
M=[J,st;-st,body.mass*eye(3)]; H=[Hm;body.S_bar];
Md=[Jd,sdt;-sdt,zeros(3)]; h=M*v+H*ed; star=crf(v);
C=Md*v+star*h; Cv=Md+star*M;
dM=zeros(6,6,n); dH=zeros(6,n,n); Ceta=zeros(6,n); Ced=Ceta;
Cv=Cv+force_right(h);
for j=1:n
    sj=get_skew_sym_matrix(body.S_bar(:,j));
    dM(:,:,j)=[inertia_from_moment(Te(:,:,j)),sj;-sj,zeros(3)];
    dMd=[inertia_from_moment(Tde(:,:,j)),zeros(3);zeros(3),zeros(3)];
    Dh=zeros(3,n);
    pairs=[3,2,2,3;1,3,3,1;2,1,1,2];
    for k=1:3
        pl=body.shapeProducts{pairs(k,1),pairs(k,2)};
        pr=body.shapeProducts{pairs(k,3),pairs(k,4)};
        Dh(k,:)=(pl(2:end,j+1)-pr(2:end,j+1)).';
    end
    dH(:,:,j)=[Dh;zeros(3,n)];
    Ceta(:,j)=dMd*v+star*(dM(:,:,j)*v+dH(:,:,j)*ed);
    Ced(:,j)=dM(:,:,j)*v+star*H(:,j);
end
for j=1:3
    unit=zeros(3,1); unit(j)=1; dW=get_skew_sym_matrix(unit);
    dW2=dW*W+W*dW; cv(:,j)=body.S_bar.'*dW*v(4:6);
    for k=1:3
        for l=1:3
            cv(:,j)=cv(:,j)+dW2(k,l)*P{k,l}+2*dW(k,l)*V{k,l};
        end
    end
end
end

