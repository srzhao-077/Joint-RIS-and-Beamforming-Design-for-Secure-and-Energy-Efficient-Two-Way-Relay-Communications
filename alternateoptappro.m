function [best_Phi,best_A] = alternateoptappro(nr,nt,Pb,sigm2,Ps,beita,itersmax_bcd,m,h1,h2,h3,g1,g2,g3,V,eta,iniphi,inia)

best_Phi = iniphi;
best_A = inia;
eta_list = 0; 
[reals, ~] = getSumrate(h1, h2, h3, g1, g2, g3, V, iniphi, inia, Ps, sigm2);
Phi = iniphi;
A = inia;
iter_bcd = 1;
thres =0.8;

    
    while iter_bcd <= itersmax_bcd
        [U1, U2, U3, W1, W2, W3, W4] = relatedWmmse(h1, h2, h3, g1, g2, g3, V, Phi, A, beita, nt, nr);
        A = updateseeA(h1, h2, h3, g1, g2, g3, V, Phi, Ps, sigm2, U1, U2, U3, W1, W2, W3, W4, Pb, nr, A, eta);
        Phi = updatephiapprocvx(m, Phi, Ps, Pb, h1, h2, h3, g1, g2, g3, V, A, beita, nr, sigm2, U1, U2, U3, W1, W2, W3, W4);
        [realsr, energy] = getSumrate(h1, h2, h3, g1, g2, g3, V, Phi, A, Ps, sigm2);
        
        if iter_bcd > 2 &&  realsr <thres * reals
            break;  
        end
                
        iter_bcd = iter_bcd + 1;
        
        if realsr/energy > eta_list
        eta_list = realsr/energy;
        best_Phi = Phi;
        best_A = A;
        end
    end
 
end