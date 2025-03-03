function [seoptsr,seopteta] = seeopt(nr,nt,Pb,sigm2,Ps,beita,m,itersmax_bcd,h1,h2,h3,g1,g2,g3,V,sroptphi,sropta)
eta=0.001;
realsr=zeros(itersmax_bcd,1);
energy=zeros(itersmax_bcd,1);
iniphi= sroptphi;
inia= sropta;
etaarray=zeros(itersmax_bcd,1); 
iter_bcd = 1;
[inisr,inie]= getSumrate(h1, h2, h3, g1, g2, g3, V, iniphi, inia, Ps, sigm2);
realsr(iter_bcd)=inisr;
energy(iter_bcd)=inie;
etaarray(iter_bcd)=inisr/inie;

while iter_bcd <= itersmax_bcd
[best_Phi,best_A] = alternateoptforsee(nr,nt,Pb,sigm2,Ps,beita,itersmax_bcd,m,h1,h2,h3,g1,g2,g3,V,eta,iniphi,inia);
[realsr(iter_bcd+1),energy(iter_bcd+1)]= getSumrate(h1, h2, h3, g1, g2, g3, V, best_Phi, best_A, Ps, sigm2);
       if energy(iter_bcd+1) > 0  
            etaarray(iter_bcd+1) = realsr(iter_bcd+1) / energy(iter_bcd+1);
       end                    
       
       if iter_bcd>1 &&  etaarray(iter_bcd+1) <etaarray(iter_bcd)
           break;
       end
        iter_bcd = iter_bcd + 1;
end

[seopteta, index] = max(etaarray);
seoptsr = realsr(index);
end