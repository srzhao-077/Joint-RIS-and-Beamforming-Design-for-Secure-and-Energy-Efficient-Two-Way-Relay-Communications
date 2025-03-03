function [ransr,raneta] = ranopt(nr,sigm2,Ps,m,h1,h2,h3,g1,g2,g3,V)
ori_xita=2*pi*randn(m,1);
Phi=diag(cos(ori_xita)+sqrt(-1)*sin(ori_xita));
A1=randi([0, 50])*(randn(nr,nr)+sqrt(-1)*(randn(nr,nr)));
[ransr,rane] = getSumrate(h1,h2,h3,g1,g2,g3,V,Phi,A1,Ps,sigm2);
raneta= ransr/rane;