function [X,dX,g]=deformation_terms(Phi,eta)
% Exact coordinate derivative of the Rodrigues transform used in Xdef.
[R,~,delta]=get_deformation_transform(Phi,eta);
theta=delta(1:3); u=delta(4:6); K=get_skew_sym_matrix(theta); s=theta.'*theta;
if s<1e-16
    aa=1-s/6+s^2/120; bb=1/2-s/24+s^2/720;
else
    t=sqrt(s); aa=sin(t)/t; bb=(1-cos(t))/s;
end
if s<1e-4
    da=-1/6+s/60-s^2/1680+s^3/90720;
    db=-1/24+s/360-s^2/13440+s^3/907200;
else
    t=sqrt(s); da=(t*cos(t)-sin(t))/(2*t^3);
    db=(t*sin(t)-2*(1-cos(t)))/(2*t^4);
end
E=R.'; U=get_skew_sym_matrix(u);
X=[E,zeros(3);-E*U,E]; dX=zeros(6,6,numel(eta)); g=zeros(6,numel(eta));
for j=1:numel(eta)
    dK=get_skew_sym_matrix(Phi(1:3,j)); ds=2*theta.'*Phi(1:3,j);
    dE=(aa*dK+bb*(dK*K+K*dK)+da*ds*K+db*ds*K*K).';
    dU=get_skew_sym_matrix(Phi(4:6,j));
    dX(:,:,j)=[dE,zeros(3);-dE*U-E*dU,dE];
    dRotation=dE*E.';
    g(:,j)=[0.5*[dRotation(3,2)-dRotation(2,3);dRotation(1,3)-dRotation(3,1);dRotation(2,1)-dRotation(1,2)]; -E*Phi(4:6,j)];
end
end

