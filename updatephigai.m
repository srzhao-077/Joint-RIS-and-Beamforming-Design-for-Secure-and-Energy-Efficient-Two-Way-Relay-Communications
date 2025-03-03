function [newPhi] = updatephigai(m, Phi, Ps, Pb, h1, h2, h3, g1, g2, g3, V, A, beita, nr, sigm2, U1, U2, U3, W1, W2, W3, W4)
    ori_xita = 2*pi*randn(m,1);
    valy_t = exp(1i * ori_xita);
    
    u_11 = g1 * h1' * A' * conj(V); 
    u_22 = g2 * h2' * A' * conj(V); 
    u_33 = g3 * h3' * A' * conj(V);
    v_11 = conj(g1) * h1.' * A * V; 
    v_22 = conj(g2) * h2.' * A * V; 
    v_33 = conj(g3) * h3.' * A * V;
    
    w = V' * A' * conj(V);
    z_12 = conj(g1) * h1.' * A * h2 * g2'; 
    z_21 = conj(g2) * h2.' * A * h1 * g1';
    z_32 = conj(g3) * h3.' * A * h2 * g2'; 
    z_31 = conj(g3) * h3.' * A * h1 * g1';
    
    Q3 = sqrt(beita) * W1 * U1' * (w.' .* (conj(g1) * g2')) + ...
         sqrt(beita) * W2 * U2' * (w.' .* (conj(g2) * g1'));
    Psi = beita * real(U1 * W1 * U1') * (w.' .* z_12) + ...
          beita * real(U2 * W2 * U2') * (w.' .* z_21) + ...
          beita * real(W4) * (w.' .* z_31) + beita * real(W4) * (w.' .* z_32) - Q3;
    
    Xi = beita * U1 * W1 * U1' * (v_11.' .* u_22) + ...
         beita * U2 * W2 * U2' * (v_22.' .* u_11) + ...
         beita * W4 * (v_33.' .* u_11) + ...
         beita * W4 * (v_33.' .* u_22);

    G1 = V * diag(g1); 
    G2 = V * diag(g2); 
    G3 = V * diag(g3);
    
    Gbar1 = conj(h1) * h1.' + m * conj(G1) * G1.'; 
    Gbar2 = conj(h2) * h2.' + m * conj(G2) * G2.'; 
    Gbar3 = conj(h3) * h3.' + m * conj(G3) * G3.';
    
    Q1 = -sqrt(beita)*W1*U1'*(G1'*conj(A)*conj(h2)+G2'*A'*conj(h1)) - ...
         sqrt(beita)*W2*U2'*(G2'*conj(A)*conj(h1)+G1'*A'*conj(h2)) - ...
         G3'*conj(A)*trace(W3)*U3 + ...
         (beita*real(U1*W1*U1')*(G2'*A'*Gbar1*A*h2 + G1'*conj(A)*Gbar2*A.'*h1) + ...
         beita*real(U2*W2*U2')*(G1'*A'*Gbar2*A*h1 + G2'*conj(A)*Gbar1*A.'*h2) + ...
         beita*real(W4)*(G1'*A'*Gbar3*A*h1 + G3'*conj(A)*Gbar1*A.'*h3) + ...
         beita*real(W4)*(G2'*A'*Gbar3*A*h2 + G3'*conj(A)*Gbar2*A.'*h3)) + ...
         (real(U1*W1*U1')*G1'*conj(A)*A.'*h1 + real(U2*W2*U2')*G2'*conj(A)*A.'*h2 + ...
         real(trace(W3))*G3'*conj(A)*(U3*U3')*A.'*h3 + real(W4)*G3'*conj(A)*A.'*h3);
    
    q2 = Ps * G1' * (A'*A) * h1 + Ps * G2' * (A'*A) * h2;
    const_fst = real(trace((Ps*Gbar1.' + Ps*Gbar2.' + sigm2*eye(nr))*(A'*A)));
    
    const_fobj = real(beita*real(U1*W1*U1')*trace(Gbar2.'*A'*Gbar1*A) + ...
                     beita*real(U2*W2*U2')*trace(Gbar1.'*A'*Gbar2*A) + ...
                     beita*real(W4)*trace(Gbar1.'*A'*Gbar3*A) + ...
                     beita*real(W4)*trace(Gbar2.'*A'*Gbar3*A) + ...
                     real(U1*W1*U1')*trace(Gbar1*(A*A')) + ...
                     real(U2*W2*U2')*trace(Gbar2*(A*A')) + ...
                     trace(W3)*trace(Gbar3*(A*(U3*U3')*A')) + ...
                     real(W4)*trace(Gbar3*(A*A'))) - ...
                     2*real(trace(sqrt(beita)*W1*U1'*h1.'*A*h2 + sqrt(beita)*W2*U2'*h2.'*A*h1 + h3.'*A*trace(W3)*conj(U3)));
    
    iter_max = 12;
    phi = zeros(m, iter_max); 
    phi(:,1) = diag(Phi);  
    
    for t = 1:iter_max-1
        [~, funcg_grad] = relatedfuncg(phi(:,t), Xi, Psi, m, Q1);
        
       H_f = (Xi'+Xi); 
        
        grad_g1 = q2;
        constraint_value = 2* grad_g1' * phi(:,t) + const_fst;
        
        A_kkt = [H_f, grad_g1; 
                grad_g1', 0];
        b = [-funcg_grad; -constraint_value];  
        
        delta = A_kkt \ b;
        delta_phi = delta(1:m);
        
        step = armijo_search(phi(:,t), funcg_grad, delta_phi, Xi, Psi, Q1);
        phi_candidate = phi(:,t) + step * delta_phi;
        phi(:,t+1) = exp(1i * angle(phi_candidate));
        
        if t > 1 && norm(phi(:,t+1) - phi(:,t)) < 1e-4
            break;
        end
    end
    
    newPhi = diag(phi(:,t+1));
    if getSumrate(h1,h2,h3,g1,g2,g3,V,newPhi,A,Ps,sigm2) < getSumrate(h1,h2,h3,g1,g2,g3,V,Phi,A,Ps,sigm2)
        newPhi = Phi; 
    end
end

function step = armijo_search(phi, grad, dir, Xi, Psi, Q1)
    alpha = 0.5;
    beta = 0.5;
    sigma = 0.1;
    max_iter = 10;
    f_curr = compute_objective(phi, Xi, Psi, Q1);
    for k = 1:max_iter
        phi_new = phi + alpha * dir;
        phi_new = exp(1i * angle(phi_new)); 
        f_new = compute_objective(phi_new, Xi, Psi, Q1);
        if f_new <= f_curr + sigma * alpha * real(grad' * dir)
            break;
        end
        alpha = alpha * beta;
    end
    step = alpha;
end

function f = compute_objective(phi, Xi, Psi, Q1)
    f = real(phi' * Xi * phi + phi' * Psi * conj(phi) + Q1' * phi);
end


function H_f = numerical_hessian(phi, func_handle, eps)
    m = length(phi);
    H_f = zeros(m, m);
    for i = 1:m
        phi_plus = phi;
        phi_plus(i) = phi_plus(i) + eps;
        [~, grad_plus] = func_handle(phi_plus);
        
        phi_minus = phi;
        phi_minus(i) = phi_minus(i) - eps;
        [~, grad_minus] = func_handle(phi_minus);
        
        H_f(:,i) = (grad_plus - grad_minus) / (2*eps);
    end
    H_f = 0.5*(H_f + H_f');
end




function step = compute_step(Xi, Psi, Q1, phi_t, dirt)
    effimuA = dirt' * Xi * dirt + dirt' * Psi * conj(dirt);
    effimuB = phi_t' * Xi * dirt + dirt' * Xi * phi_t ...
             + phi_t' * Psi * conj(dirt) + dirt' * Psi * conj(phi_t) ...
             + dirt' * Q1;
    effimuC = phi_t' * Xi * phi_t + phi_t' * Psi * conj(phi_t) + phi_t' * Q1;

   
    step = bisection_method(effimuA, effimuB, effimuC);
end




function step = bisection_method(effimuA, effimuB, effimuC)
    alpha = 0.618;
    a = 0; b = 1;
    accu = 10^(-4);
    shrinkDirection = 3;

    while b - a > accu
        if shrinkDirection == 3
            lambda_k = a + (1 - alpha) * (b - a);
            mu_k = a + alpha * (b - a);
        else
            lambda_k = mu_k;
            mu_k = a + alpha * (b - a);
        end

        orig_lambda_k = 2 * real(effimuA * lambda_k^2 + effimuB * lambda_k + effimuC);
        orig_mu_k = 2 * real(effimuA * mu_k^2 + effimuB * mu_k + effimuC);

        if abs(orig_lambda_k - orig_mu_k) <= accu
            break;
        end
        if orig_lambda_k <= orig_mu_k
            b = mu_k;
            shrinkDirection = 1;
        else
            a = lambda_k;
            shrinkDirection = 0;
        end
    end
    step = (b + a) / 2;
end
