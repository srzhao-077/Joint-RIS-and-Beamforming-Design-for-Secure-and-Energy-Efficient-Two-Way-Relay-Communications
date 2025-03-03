 
function [h1,h2,h3,g1,g2,g3,V] = randchannel(nr,nt,m)

K=10;
x_bs=0; y_bs=0;
x_ris=40; y_ris=20;
x_s1=60; y_s1=10.;
x_s2=60; y_s2=-10;
x_s3=80; y_s3=0;



d_h1=sqrt((x_bs-x_s1)^2+(y_bs-y_s1)^2);
d_h2=sqrt((x_bs-x_s2)^2+(y_bs-y_s2)^2);
d_h3=sqrt((x_bs-x_s3)^2+(y_bs-y_s3)^2);
d_g1=sqrt((x_ris-x_s1)^2+(y_ris-y_s1)^2);
d_g2=sqrt((x_ris-x_s2)^2+(y_ris-y_s2)^2);
d_g3=sqrt((x_ris-x_s3)^2+(y_ris-y_s3)^2);
d_v=sqrt((x_bs-x_ris)^2+(y_bs-y_ris)^2);

h1=sqrt(getPathloss(3.5,d_h1))*(sqrt(1/2)*(randn(nr,nt)+sqrt(-1)*randn(nr,nt)));
h2=sqrt(getPathloss(3.5,d_h2))*(sqrt(1/2)*(randn(nr,nt)+sqrt(-1)*randn(nr,nt)));
h3=sqrt(getPathloss(3.5,d_h3))*(sqrt(1/2)*(randn(nr,nt)+sqrt(-1)*randn(nr,nt)));
aR=getRandaR(m); 
g1=sqrt(getPathloss(2.5,d_g1))*(sqrt(K/(K+1))*aR'+sqrt(1/(K+1))*sqrt(1/2)*(randn(m,nt)+sqrt(-1)*randn(m,nt)));
aR=getRandaR(m);
g2=sqrt(getPathloss(2.5,d_g2))*(sqrt(K/(K+1))*aR'+sqrt(1/(K+1))*sqrt(1/2)*(randn(m,nt)+sqrt(-1)*randn(m,nt)));
aR=getRandaR(m);
g3=sqrt(getPathloss(2.5,d_g3))*(sqrt(K/(K+1))*aR'+sqrt(1/(K+1))*sqrt(1/2)*(randn(m,nt)+sqrt(-1)*randn(m,nt)));
aB=getRandaB(nr);
V=sqrt(getPathloss(2.2,d_v))*(sqrt(K/(K+1))*aB'+sqrt(1/(K+1))*sqrt(1/2)*(randn(nr,m)+sqrt(-1)*randn(nr,m)));



end



function loss = getPathloss(erfa,d)

if d<=1
        lossdb=0;
    else
        if erfa<=2.2
            lossdb=0;
        else 
          lossdb=-30-10*erfa*log10(d);
        end
    end

    loss=10^(lossdb/10);



    
end



function [aB] = getRandaB(n)

    xita=pi*unifrnd(-1,1);
    for i=1:n
        aB(i)=exp((i-1)*(-1)*pi*sqrt(-1)*sin(xita));
    end
end

function [aR] = getRandaR(m)
    xita=pi*unifrnd(-1,1); 
    psi=25/180*pi*unifrnd(-1,1);
    mh=m/2; 
    mv=2; 
    
    for i=1:mh
        aRh(i)=exp((i-1)*(-1)*pi*sqrt(-1)*cos(psi)*sin(xita));
    end
    for i=1:mv
        aRv(i)=exp((i-1)*pi*sqrt(-1)*cos(psi)*cos(xita));
    end
    aR=kron(aRv,aRh);
end






