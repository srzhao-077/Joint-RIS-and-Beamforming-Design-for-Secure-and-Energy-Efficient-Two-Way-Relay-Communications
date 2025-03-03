
function [ranphioptsr,ranphiopteta] = ranphiopt(nr,nt,Pb,sigm2,Ps,beita,m,itersmax_bcd,h1,h2,h3,g1,g2,g3,V)
realsr=zeros(itersmax_bcd,1);
energy=zeros(itersmax_bcd,1);
eta=0.001;
etaarray = eta; 
iter_bcd = 1;
inia= randi([0, 100])*(randn(nr,nr)+sqrt(-1)*(randn(nr,nr)));
ori_xita=2*pi*randn(m,1);
best_Phi=diag(cos(ori_xita)+sqrt(-1)*sin(ori_xita));
iniphi= best_Phi;
while iter_bcd <= itersmax_bcd
best_A = ranphaseee(nr,nt,Pb,sigm2,Ps,beita,itersmax_bcd,m,h1,h2,h3,g1,g2,g3,V,eta,iniphi,inia);
[realsr(iter_bcd+1),energy(iter_bcd+1)]= getSumrate(h1, h2, h3, g1, g2, g3, V, best_Phi, best_A, Ps, sigm2);
       if energy(iter_bcd+1) > 0  
            eta = realsr(iter_bcd+1) / energy(iter_bcd+1);
       end
       etaarray=[etaarray,eta]; 
       if iter_bcd>1 &&  etaarray(iter_bcd+1) <etaarray(iter_bcd)
           break;
       end
        iter_bcd = iter_bcd + 1;
end

ranphioptsr =  max(realsr);
ranphiopteta=  max(etaarray);
end