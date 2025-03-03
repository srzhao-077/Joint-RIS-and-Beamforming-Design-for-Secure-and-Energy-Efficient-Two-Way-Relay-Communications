function [best_A] = ranphaseee(nr,nt,Pb,sigm2,Ps,beita,itersmax_bcd,m,h1,h2,h3,g1,g2,g3,V,eta,iniphi,inia)


best_A = 0;
num_groups = 3;
realsr_list = zeros(num_groups, 1); 
max_realsr = 0 ;


lhs_A = latin_hypercube_sampling(nr, nr, num_groups);


for group = 1:num_groups

    Phi = iniphi;
    A = lhs_A(:, group);  
    A = inia + A * 0.05;  
    
    [Phi, A] = apply_initial_disturbance(Phi, A); 
    iter_bcd = 1;
    while iter_bcd <= itersmax_bcd
        [U1, U2, U3, W1, W2, W3, W4] = relatedWmmse(h1, h2, h3, g1, g2, g3, V, Phi, A, beita, nt, nr);
        A = updateAee(h1, h2, h3, g1, g2, g3, V, Phi, Ps, sigm2, U1, U2, U3, W1, W2, W3, W4, Pb, nr, A, eta);
        [realsr, energy] = getSumrate(h1, h2, h3, g1, g2, g3, V, Phi, A, Ps, sigm2);
        
        if iter_bcd > 2 && abs(realsr-eta*energy - realsr_list(group)) < 1e-3
            break;  
        end
        
        if iter_bcd > 2 && realsr < 0.5
            break;
        end
        
        iter_bcd = iter_bcd + 1;
    end
    realsr_list(group) = realsr-eta*energy;
  
    if realsr-eta*energy > max_realsr
        max_realsr = realsr-eta*energy;
        best_A = A;
    end
end

end

function [Phi, A] = apply_initial_disturbance(Phi, A)
   
    delta = 0.01;  
    
    Phi = Phi + delta * randn(size(Phi)); 
    A = A + delta * randn(size(A));  
end

function [U1,U2,U3,W1,W2,W3,W4] = relatedWmmse(h1,h2,h3,g1,g2,g3,V,Phi,A,beita,nt,nr)



hg1=h1+V*Phi*g1; hg2=h2+V*Phi*g2; hg3=h3+V*Phi*g3;
H1B=hg1.'*A; H2B=hg2.'*A; H3B=hg3.'*A;% 
H12=hg1.'*A*hg2; H21=hg2.'*A*hg1; H31=hg3.'*A*hg1; H32=hg3.'*A*hg2; 

% Calculate l1, l2, l3
U1 = sqrt(beita) * (inv( hg1.' * A * (beita * (hg2 * hg2') + eye(nr)) * A' * conj(hg1) + 1 ))' * hg1.' * A * hg2;
U2 = sqrt(beita) * (inv( hg2.' * A * (beita * (hg1 * hg1') + eye(nr)) * A' * conj(hg2) + 1 ))' * hg2.' * A * hg1;
U3 = inv(eye(nr) + A.' * (hg3*hg3') * conj(A)) * A.' * hg3;




% Calculate f1, f2, f3, f4
W1 = 1 + beita * (H12'* H12)/(H1B*H1B'+1) ;
W2 = 1 + beita * (H21'* H21)/(H2B*H2B'+1) ;
W3 = eye(nr) + H3B'*H3B;
W4 =1 / (1 + H3B * H3B' + beita * (H31*H31'+H32*H32'));

end


function [newA] = updateAee(h1,h2,h3,g1,g2,g3,V,Phi,Ps,sigm2,U1,U2,U3,W1,W2,W3,W4,Pb,nr,oldA,eta)


hg1=h1+V*Phi*g1; hg2=h2+V*Phi*g2; hg3=h3+V*Phi*g3;
beita=Ps/sigm2;
linyong = eta.*(Ps.*hg1*hg1'+Ps.*hg2*hg2'+sigm2*eye(nr));
T1=kron((beita*(hg2*hg2')+eye(nr)).',(conj(hg1)*real(U1*W1*U1')*hg1.')) ...
    + kron((beita*(hg1*hg1')+eye(nr)).',(conj(hg2)*real(U2*W2*U2')*hg2.'))...
    + kron((conj(U3)*U3.').',trace(W3)*(conj(hg3)*hg3.'))...
    + kron((beita*(hg1*hg1'+hg2*hg2')+eye(nr)).',(conj(hg3)*W4*hg3.')) + kron(eye(nr),linyong);
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



end


function lhs = latin_hypercube_sampling(dim, size, num_samples)
    lhs = zeros(dim, num_samples);
    for i = 1:num_samples
        lhs(:, i) = (lhsdesign(dim, 1) - 0.5) * size;
    end
end
