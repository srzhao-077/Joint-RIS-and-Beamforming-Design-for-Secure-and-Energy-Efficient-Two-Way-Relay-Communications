clear;
nr = 3;
nt = 1; 
Pb = 10^(10/10); 
sigm2 = 10^(-104/10); 
Ps = 10^(0/10); 
beita = Ps/sigm2;
num_iterations = 10000;
itersmax_bcd = 5;

m_values = [10, 20, 30, 40];  
num_m = length(m_values);

total_time1 = zeros(1, num_m); 
total_time2 = zeros(1, num_m);

for idx = 1:num_m
    m = m_values(idx);
    
    seoptsubsr_total = 0;
    seoptsubeta_total = 0;
    seoptsr_total = 0;
    seopteta_total = 0;
    
    for iter = 1:num_iterations 
        [h1,h2,h3,g1,g2,g3,V] = randchannel(nr,nt,m);
        [initializephi,initializea] = initialize(nr,nt,Pb,sigm2,Ps,beita,m,itersmax_bcd,h1,h2,h3,g1,g2,g3,V);
        
        tic;
        [seoptsubsr, seoptsubeta] = seeoptsub(nr,nt,Pb,sigm2,Ps,beita,m,itersmax_bcd,h1,h2,h3,g1,g2,g3,V,initializephi,initializea);
        
        seoptsubsr_total = seoptsubsr_total + seoptsubsr;
        seoptsubeta_total = seoptsubeta_total + seoptsubeta;
        elapsedTime1 = toc; 
        total_time1(idx) = total_time1(idx) + elapsedTime1;
        
        tic;
        [seoptsr, seopteta] = seeoptappro(nr,nt,Pb,sigm2,Ps,beita,m,itersmax_bcd,h1,h2,h3,g1,g2,g3,V,initializephi,initializea);
        
        seoptsr_total = seoptsr_total + seoptsr;
        seopteta_total = seopteta_total + seopteta;
        elapsedTime2 = toc;
        total_time2(idx) = total_time2(idx) + elapsedTime2;
    end
    
    
    
end
total_time1 = total_time1/num_iterations;
total_time2 = total_time2/num_iterations;
