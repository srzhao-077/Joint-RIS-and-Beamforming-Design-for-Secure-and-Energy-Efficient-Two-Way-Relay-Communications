
function [const_g,funcgd] = relatedfuncg(phit,Xi,Psi,m,vcoeffi1)
vcoeffi2=(Psi.'+Psi)*conj(phit);
vcoeffi3=(Xi'+Xi)*phit;
const_g=real(-phit'*Psi.'*conj(phit))+real(-phit.'*Xi.'*conj(phit));
funcgd=(vcoeffi1+vcoeffi2+vcoeffi3);
end