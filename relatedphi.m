function [step,rate] = relatedphi(Xi,Psi,Q1,phit,dirt)
effimuA=dirt'*Xi*dirt+dirt'*Psi*conj(dirt);
effimuB=phit'*Xi*dirt+dirt'*Xi*phit...
         +phit'*Psi*conj(dirt)+dirt'*Psi*conj(phit)...
         +dirt'*Q1;
effimuCsty=phit'*Xi*phit+phit'*Psi*conj(phit)+phit'*Q1;
accu=10^(-4);
count=0;
alpha=0.618;
a=0; b=1;
shrinkDirection=3;

while b-a >accu 

    if shrinkDirection==3
        lambda_k = a + (1 - alpha) * (b - a);
        mu_k = a + alpha * (b - a);
    else
        if shrinkDirection==1
            mu_k=lambda_k;
            lambda_k = a + (1 - alpha) * (b - a);
        else
            if shrinkDirection==0
                lambda_k=mu_k;
                mu_k = a + alpha * (b - a);
            end
        end
    end

    orig_lambda_k=2*real(effimuA*lambda_k^2+effimuB*lambda_k+effimuCsty);
    orig_mu_k=2*real(effimuA*mu_k^2+effimuB*mu_k+effimuCsty);

    if  abs(orig_lambda_k-orig_mu_k)<=accu
        phit1=phit;
        break;
    end
    if(orig_lambda_k)<=(orig_mu_k)
        b=mu_k;
        shrinkDirection=1;
    else
        a=lambda_k;
        shrinkDirection=0;
    end
    count = count+1;
    step=(b+a)/2;
    phit1=phit+step*dirt;
end
step=(b+a)/2;
rate=2*real(effimuA*step^2+effimuB*step+effimuCsty);

end