function [newA] = updateA(h1,h2,h3,g1,g2,g3,V,Phi,Ps,sigm2,U1,U2,U3,W1,W2,W3,W4,Pb,nr,oldA)

hg1=h1+V*Phi*g1; hg2=h2+V*Phi*g2; hg3=h3+V*Phi*g3;
beita=Ps/sigm2;

T1=kron((beita*(hg2*hg2')+eye(nr)).',(conj(hg1)*real(U1*W1*U1')*hg1.')) ...
    + kron((beita*(hg1*hg1')+eye(nr)).',(conj(hg2)*real(U2*W2*U2')*hg2.'))...
    + kron((conj(U3)*U3.').',trace(W3)*(conj(hg3)*hg3.'))...
    + kron((beita*(hg1*hg1'+hg2*hg2')+eye(nr)).',(conj(hg3)*W4*hg3.'));
X1=sqrt(beita)*conj(hg1)*U1*W1'*hg2'+sqrt(beita)*conj(hg2)*U2*W2'*hg1'+conj(hg3)*trace(W3)*U3.';

T2=Ps*kron((hg1*hg1').',eye(nr))+Ps*kron((hg2*hg2').',eye(nr))+kron(sigm2*eye(nr),eye(nr));
vk =vec(X1);
T1_inv = eye(nr^2)/T1;
T2_inv = eye(nr^2)/T2;

co_a=real(vk' * (T1_inv' * T2 * T1_inv) *vk)-Pb;
co_b=2*real(vk' * T2_inv' * T2 * T1_inv * vk);
co_c=real(vk' * T2_inv' * T2 * T2_inv * vk);
co_t=(co_b)^2-4*(co_a)*(co_c);
if co_t<0
    lambda=0;
    a1=inv(T1)*vk;
    newA=oldA;
else
    lambda=(-co_b-sqrt(co_t))/(2*co_a);
    a2=inv(T1+lambda*T2)*vk;
    newA=reshape(a2,[nr,nr]);
end


Val1=getSumrate(h1,h2,h3,g1,g2,g3,V,Phi,oldA,Ps,sigm2);
Val2=getSumrate(h1,h2,h3,g1,g2,g3,V,Phi,newA,Ps,sigm2);
if Val1>Val2
    newA=oldA;
end
end