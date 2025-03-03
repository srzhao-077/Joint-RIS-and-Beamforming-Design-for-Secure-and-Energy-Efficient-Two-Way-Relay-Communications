function [sumrate,energy] = getSumrate(h1,h2,h3,g1,g2,g3,V,Phi,A,Ps,sigm2)

hg1=h1+V*Phi*g1; hg2=h2+V*Phi*g2; hg3=h3+V*Phi*g3;
H1B=hg1.'*A; H2B=hg2.'*A; H3B=hg3.'*A;
H12=hg1.'*A*hg2; H21=hg2.'*A*hg1; 
H31=hg3.'*A*hg1; H32=hg3.'*A*hg2; 



R1R2=1/2*log2(real(1+Ps*(H12*H12')/(sigm2+sigm2*H1B*H1B'))) ...
    +1/2*log2(real(1+Ps*(H21*H21')/(sigm2+sigm2*H2B*H2B')));
R3=1/2*log2(real(1+Ps*(H31*H31'+H32*H32')/(sigm2+sigm2*H3B*H3B')));



sumrate=max(0,R1R2-R3);
energy = 2 * Ps + Ps * norm(A * hg1)^2 +  Ps * norm(A * hg2)^2 + sigm2 * norm(A)^2;
end

