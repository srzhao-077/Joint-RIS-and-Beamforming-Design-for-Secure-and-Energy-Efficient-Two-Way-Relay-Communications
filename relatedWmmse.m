function [U1,U2,U3,W1,W2,W3,W4] = relatedWmmse(h1,h2,h3,g1,g2,g3,V,Phi,A,beita,nt,nr)

hg1=h1+V*Phi*g1; hg2=h2+V*Phi*g2; hg3=h3+V*Phi*g3;
H1B=hg1.'*A; H2B=hg2.'*A; H3B=hg3.'*A;
H12=hg1.'*A*hg2; H21=hg2.'*A*hg1; H31=hg3.'*A*hg1; H32=hg3.'*A*hg2; 



% Calculate l1, l2, l3
U1 = sqrt(beita) * (inv( hg1.' * A * (beita * (hg2 * hg2') + eye(nr)) * A' * conj(hg1) + 1 ))' * hg1.' * A * hg2;
U2 = sqrt(beita) * (inv( hg2.' * A * (beita * (hg1 * hg1') + eye(nr)) * A' * conj(hg2) + 1 ))' * hg2.' * A * hg1;
U3 = inv(eye(nr) + A.' * (hg3*hg3') * conj(A)) * A.' * hg3;




% Calculate f1, f2, f3, f4
W1 = 1 + beita * (H12'* H12)/(H1B*H1B'+1) ;
W2 = 1 + beita * (H21'* H21)/(H2B*H2B'+1) ;
W3 = eye(nr) + H3B'*H3B;
W4 =1 / (1 + H3B * H3B' + beita * (H31*H31'+H32*H32'));

end