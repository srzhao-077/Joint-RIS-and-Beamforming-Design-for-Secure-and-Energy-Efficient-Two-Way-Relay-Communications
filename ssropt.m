
function [sroptsr,sropteta,sroptphi,sropta] = ssropt(nr,nt,Pb,sigm2,Ps,beita,m,itersmax_bcd,h1,h2,h3,g1,g2,g3,V)

[best_Phi,best_A] = alternateopt(nr,nt,Pb,sigm2,Ps,beita,itersmax_bcd,m,h1,h2,h3,g1,g2,g3,V);
[realsr,energy]= getSumrate(h1, h2, h3, g1, g2, g3, V, best_Phi, best_A, Ps, sigm2);
        
sroptsr =  realsr;
sropteta=  realsr/energy;
sroptphi= best_Phi;
sropta  = best_A;
end