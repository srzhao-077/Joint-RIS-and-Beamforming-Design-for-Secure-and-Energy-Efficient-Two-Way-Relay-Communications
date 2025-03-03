function [newPhi] = updatephiee_grad(m, Phi, Ps, Pb, h1, h2, h3, g1, g2, g3, V, A, beita, nr, sigm2, U1, U2, U3, W1, W2, W3, W4, eta)

Oldphi =Phi;

eta0 = eta;
u_11 = g1 * h1' * A' * conj(V); u_22 = g2 * h2' * A' * conj(V); u_33 = g3 * h3' * A' * conj(V);
v_11 = conj(g1) * h1.' * A * V; v_22 = conj(g2) * h2.' * A * V; v_33 = conj(g3) * h3.' * A * V;
w = V' * A' * conj(V);
z_12 = conj(g1) * h1.' * A* h2 * g2'; z_21 = conj(g2) * h2.' * A* h1 * g1';
z_32 = conj(g3) * h3.' * A* h2 * g2'; z_31 = conj(g3) * h3.' * A* h1 * g1';
Q3 = sqrt(beita)* W1 * U1' *( w.' .* (conj(g1) * g2') ) ...
    + sqrt(beita)* W2 * U2' *( w.' .* (conj(g2) * g1') ) ;
Psi = beita * real(U1 * W1 * U1') *( w.' .* z_12 )...
    + beita * real(U2 * W2 * U2') *( w.' .* z_21 )...
    + beita * real(W4) *( w.' .* z_31 )...
    + beita * real(W4) *( w.' .* z_32 ) - Q3;
Xi = beita * U1 * W1 * U1' *( v_11.' .* u_22 )...
    + beita * U2 * W2 * U2' *( v_22.' .* u_11 )...
    + beita * W4 * ( v_33.' .* u_11 )...
    + beita * W4 * ( v_33.' .* u_22 );
G1=V*diag(g1); G2=V*diag(g2); G3=V*diag(g3);
Gbar1=conj(h1)*h1.'+m*conj(G1)*G1.'; Gbar2=conj(h2)*h2.'+m*conj(G2)*G2.'; Gbar3=conj(h3)*h3.'+m*conj(G3)*G3.';
Q1=-sqrt(beita)*W1*U1'*(G1'*conj(A)*conj(h2) + G2'*A'*conj(h1))...
    -sqrt(beita)*W2*U2'*(G2'*conj(A)*conj(h1) + G1'*A'*conj(h2))...
    -G3'*conj(A)*trace(W3)*U3...
    +(beita*real(U1 * W1 * U1')*(G2'*A'*Gbar1*A*h2+G1'*conj(A)*Gbar2*A.'*h1)...
    +beita*real(U2 * W2 * U2')*(G1'*A'*Gbar2*A*h1+G2'*conj(A)*Gbar1*A.'*h2)...
    +beita*real(W4)*(G1'*A'*Gbar3*A*h1+G3'*conj(A)*Gbar1*A.'*h3)...
    +beita*real(W4)*(G2'*A'*Gbar3*A*h2+G3'*conj(A)*Gbar2*A.'*h3))...
    +(real(U1 * W1 * U1')*G1'*conj(A)*A.'*h1+real(U2 * W2 * U2')*G2'*conj(A)*A.'*h2...
    +real(trace(W3))*G3'*conj(A)*(U3*U3')*A.'*h3+real(W4)*G3'*conj(A)*A.'*h3);
q2=Ps*G1'*(A'*A)*h1+Ps*G2'*(A'*A)*h2;
const_fst=real(trace((Ps*Gbar1.'+Ps*Gbar2.'+sigm2*eye(nr))*(A'*A)));
const_fobj=real(beita*real(U1*W1*U1')*trace(Gbar2.'*A'*Gbar1*A) ...
    +beita*real(U2 * W2 * U2')*trace(Gbar1.'*A'*Gbar2*A)...
    +beita*real(W4)*trace(Gbar1.'*A'*Gbar3*A)+beita*real(W4)*trace(Gbar2.'*A'*Gbar3*A)...
    +real(U1*W1*U1')*trace(Gbar1*(A*A'))+real(U2 * W2 * U2')*trace(Gbar2*(A*A'))...
    +trace(W3)*trace(Gbar3*(A*(U3*U3')*A'))+real(W4)*trace(Gbar3*(A*A')))...
    -2*real(trace(sqrt(beita)*W1*U1'*h1.'*A*h2+sqrt(beita)*W2*U2'*h2.'*A*h1+h3.'*A*trace(W3)*conj(U3)));

    iter_max = 10;
    phi = zeros(m, iter_max);
    phi_y = zeros(m, iter_max + 1);
    temp2sumrate = zeros(iter_max, 1);
    
    for k = 1:m
        phi(k, 1) = Phi(k, k);
    end
    phi_y(:, 1) = phi(:, 1);
    
    t = 1;
    while (t < iter_max)
        [const_g, funcg_grad] = relatedfuncg(phi(:, t), Xi, Psi, m, Q1); 
        

        grad_f = 2 * real(funcg_grad) + const_g; 
        
   
        grad_step = -eta0 * grad_f;  
        valy_t = phi(:, t) + grad_step;  
        

        for k = 1:m
            valy_t(k) = valy_t(k) / norm(valy_t(k)); 
        end
        
        cc = 1;
        while 2 * real(valy_t' * q2) + const_fst > Pb && cc < 100

            valy_t = valy_t - grad_step * 0.5;  
            
            cc=cc+1;
        end
        

        phi_y(:, t + 1) = valy_t;
        

        dirt = phi_y(:, t + 1) - phi(:, t);
        

        [step, temp2sumrate(t + 1)] = relatedphi(Xi, Psi, Q1, phi(:, t), dirt); 
        phi(:, t + 1) = phi(:, t) + step * dirt;
        
        t = t + 1;
        
        if abs((temp2sumrate(t) - temp2sumrate(t - 1)) / temp2sumrate(t - 1)) <= 10^(-4)
            break;
        end
    end
    
    phit = phi(:, t);
    newPhi = diag(phit);
    
[Val1,Va33]=getSumrate(h1,h2,h3,g1,g2,g3,V,Oldphi,A,Ps,sigm2);
[Val2,Va44]=getSumrate(h1,h2,h3,g1,g2,g3,V,newPhi,A,Ps,sigm2);
if Val1/Va33>Val2/Va44
    newPhi=Phi;
end

end