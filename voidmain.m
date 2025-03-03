function [seopteta_avg,sropteta_avg,seoptsubeta_avg,ranphiopteta_avg,raneta_avg,cceta_avg,seoptsr_avg,sroptsr_avg,seoptsubsr_avg,ranphioptsr_avg,ransr_avg,ccsr_avg] = voidmain(nr,nt,Pb,sigm2,Ps,beita,m,itersmax_bcd)

num_iterations = 10000;

ccsr_total = 0;
cceta_total = 0;
ransr_total = 0;
raneta_total = 0;
ranphioptsr_total = 0;
ranphiopteta_total = 0;
sroptsr_total = 0;
sropteta_total = 0;
seoptsubsr_total = 0;
seoptsubeta_total = 0;
seoptsr_total = 0;
seopteta_total = 0;

for iter = 1:num_iterations 
   
[h1,h2,h3,g1,g2,g3,V] = randchannel(nr,nt,m);

[ccsr,cceta] = ccopt(nr,sigm2,Ps,m,h1,h2,h3,g1,g2,g3,V) ;
ccsr_total = ccsr_total + ccsr;
cceta_total = cceta_total + cceta;

[ransr,raneta] = ranopt(nr,sigm2,Ps,m,h1,h2,h3,g1,g2,g3,V);
ransr_total = ransr_total + ransr;
raneta_total = raneta_total + raneta;

[ranphioptsr,ranphiopteta] = ranphiopt(nr,nt,Pb,sigm2,Ps,beita,m,itersmax_bcd,h1,h2,h3,g1,g2,g3,V);
ranphioptsr_total = ranphioptsr_total + ranphioptsr;
ranphiopteta_total = ranphiopteta_total + ranphiopteta;
 
[sroptsr,sropteta,sroptphi,sropta] = ssropt(nr,nt,Pb,sigm2,Ps,beita,m,itersmax_bcd,h1,h2,h3,g1,g2,g3,V);
sroptsr_total = sroptsr_total + sroptsr;
sropteta_total = sropteta_total + sropteta;

[seoptsubsr,seoptsubeta] = seeoptsub(nr,nt,Pb,sigm2,Ps,beita,m,itersmax_bcd,h1,h2,h3,g1,g2,g3,V,sroptphi,sropta);
seoptsubsr_total = seoptsubsr_total + seoptsubsr;
seoptsubeta_total = seoptsubeta_total + seoptsubeta;

[seoptsr,seopteta] = seeopt(nr,nt,Pb,sigm2,Ps,beita,m,itersmax_bcd,h1,h2,h3,g1,g2,g3,V,sroptphi,sropta);
seoptsr_total = seoptsr_total + seoptsr;
seopteta_total = seopteta_total + seopteta;
end



ccsr_avg = ccsr_total / num_iterations * 10;
cceta_avg = cceta_total / num_iterations * 10;
ransr_avg = ransr_total / num_iterations * 10;
raneta_avg = raneta_total / num_iterations * 10;
ranphioptsr_avg = ranphioptsr_total / num_iterations * 10;
ranphiopteta_avg = ranphiopteta_total / num_iterations * 10;
sroptsr_avg = sroptsr_total / num_iterations * 10;
sropteta_avg = sropteta_total / num_iterations * 10;
seoptsubsr_avg = seoptsubsr_total / num_iterations * 10;
seoptsubeta_avg = seoptsubeta_total / num_iterations * 10;
seoptsr_avg = seoptsr_total / num_iterations * 10;
seopteta_avg = seopteta_total / num_iterations * 10;


