function [Mrr,Mrf,Cv,Ce,cv,ce] = velocity_terms(body,v,eta,eta_d)
% Contract each cached shape product once for both derivative families.
n = body.nEta; qa = [1;eta]; qad = [0;eta_d];
T = zeros(3); Td = T; Te = zeros(3,3,n); P = cell(3); V = cell(3);
omega = v(1:3); u = v(4:6); W = get_skew_sym_matrix(omega);
cv = zeros(n,6); ce = body.D_ff;
cv(:,4:6) = body.S_bar.'*W;
for k = 1:3
    for l = 1:3
        product = body.shapeProducts{k,l}; pq = product*qa; pqd = product*qad;
        T(k,l) = qa.'*pq;
        Td(k,l) = qad.'*pq+qa.'*pqd;
        Te(k,l,:) = reshape(pq(2:end)+product(:,2:end).'*qa,1,1,n);
        P{k,l} = reshape(pq(2:end),n,1); V{k,l} = reshape(pqd(2:end),n,1);
        ce = ce+2*W(k,l)*product(2:end,2:end);
    end
end
J = trace(T)*eye(3)-T; 
Jd = trace(Td)*eye(3)-Td; 
Jd = (Jd+Jd.')/2;
H = [(P{3,2}-P{2,3}).'; 
    (P{1,3}-P{3,1}).'; 
    (P{2,1}-P{1,2}).'];
s = body.S_t0+body.S_bar*eta; 
sd = body.S_bar*eta_d;
st = get_skew_sym_matrix(s); 
sdt = get_skew_sym_matrix(sd);
Mrr = [J,st;-st,body.mass*eye(3)]; 
Mrf = [H;body.S_bar];
Md = [Jd,sdt;-sdt,zeros(3)]; 
h = Mrr*v+Mrf*eta_d;
Vstar = crf(v); 
Cv = Md+Vstar*Mrr; 
Ce = zeros(6,n);
for j = 1:6
    unit = zeros(6,1); unit(j) = 1;
    Cv(:,j) = Cv(:,j)+crf(unit)*h;
end
for j = 1:n
    Jj = trace(Te(:,:,j))*eye(3)-Te(:,:,j); Jj = (Jj+Jj.')/2;
    sj = get_skew_sym_matrix(body.S_bar(:,j));
    Ce(:,j) = [Jj,sj;-sj,zeros(3)]*v+Vstar*Mrf(:,j);
end
for j = 1:3
    unit = zeros(3,1); 
    unit(j) = 1; 
    dW = get_skew_sym_matrix(unit);
    dW2 = dW*W+W*dW;
    cv(:,j) = body.S_bar.'*dW*u;
    for k = 1:3
        for l = 1:3
            cv(:,j) = cv(:,j)+dW2(k,l)*P{k,l}+2*dW(k,l)*V{k,l};
        end
    end
end
end