function [Phi, A] = apply_initial_disturbance(Phi, A)
    delta = 0.01;  
    
    Phi = Phi + delta * randn(size(Phi)); 

    A = A + delta * randn(size(A)); 
end