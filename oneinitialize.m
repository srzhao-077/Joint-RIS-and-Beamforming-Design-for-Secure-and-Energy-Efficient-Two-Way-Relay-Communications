function [best_Phi,best_A] = oneinitialize(nr,nt,Pb,sigm2,Ps,beita,itersmax_bcd,m,h1,h2,h3,g1,g2,g3,V)

best_Phi = 0;
best_A = 0;
num_groups = 1;
realsr_list = zeros(num_groups, 1); 
max_realsr = 0 ;

lhs_A = latin_hypercube_sampling(nr, nr, num_groups);
for group = 1:num_groups
    
    Phi = eye(m);  
    A = lhs_A(:, group);
    A = 10*(randn(nr,nr)+sqrt(-1)*(randn(nr,nr))) + A;  
     
    [Phi, A] = apply_initial_disturbance(Phi, A);  
    iter_bcd = 1;
    while iter_bcd <= itersmax_bcd
        [U1, U2, U3, W1, W2, W3, W4] = relatedWmmse(h1, h2, h3, g1, g2, g3, V, Phi, A, beita, nt, nr);
        A = updateA(h1, h2, h3, g1, g2, g3, V, Phi, Ps, sigm2, U1, U2, U3, W1, W2, W3, W4, Pb, nr, A);
        Phi = updatephi(m, Phi, Ps, Pb, h1, h2, h3, g1, g2, g3, V, A, beita, nr, sigm2, U1, U2, U3, W1, W2, W3, W4);
        [realsr, ~] = getSumrate(h1, h2, h3, g1, g2, g3, V, Phi, A, Ps, sigm2);
        
        if iter_bcd > 2 && abs(realsr - realsr_list(group)) < 1e-3
            break;  
        end
        
        if iter_bcd > 2 && realsr < 0.1
            break;  
        end
        
        iter_bcd = iter_bcd + 1;
    end
    realsr_list(group) = realsr;

    
    if realsr > max_realsr
        max_realsr = realsr;
        best_Phi = Phi;
        best_A = A;
    end
end

end
