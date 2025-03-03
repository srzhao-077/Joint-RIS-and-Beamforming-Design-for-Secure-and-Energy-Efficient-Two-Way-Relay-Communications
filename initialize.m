
function [initializephi,initializea] = initialize(nr,nt,Pb,sigm2,Ps,beita,m,itersmax_bcd,h1,h2,h3,g1,g2,g3,V)

[best_Phi,best_A] = oneinitialize(nr,nt,Pb,sigm2,Ps,beita,itersmax_bcd,m,h1,h2,h3,g1,g2,g3,V);
initializephi= best_Phi;
initializea  = best_A;
end