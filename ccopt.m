function [ccsr,cceta] = ccopt(nr,sigm2,Ps,m,h1,h2,h3,g1,g2,g3,V)    
A2=eye(nr);
Phi2=zeros(m);
[ccsr,ccne] = getSumrate(h1,h2,h3,g1,g2,g3,V,Phi2,A2,Ps,sigm2);
cceta= ccsr/ccne;