# Joint-RIS-and-Beamforming-Design-for-Secure-and-Energy-Efficient-Two-Way-Relay-Communications
This is the code package related to the following scientific article: "Joint RIS and Beamforming Design for Secure and Energy-Efficient Two-Way Relay Communications,"  Submitted to IEEE Transactions on Mobile Computing. The package contains a simulation environment, based on Matlab.

# Abstract of Article
This paper examines the enhancement of secrecy energy efficiency (SEE) in a reconfigurable intelligent surface (RIS)-assisted two-way relay (TWR) system. We first establish a theoretical model for the system’s secrecy rate, energy consumption, and SEE, and formulate the SEE maximization problem through the joint design of the RIS phase shifts and beamforming matrix. Using techniques such as weighted minimum mean square error (WMMSE), alternating optimization, and the augmented Lagrange method, we then develop a theoretical framework that identifies locally optimal solutions for the RIS and beamforming settings under unit-modulus and power constraints. The proposed framework is also shown to be applicable to solving the system’s secrecy rate maximization
problem. To address the computational complexity involved in optimizing the RIS phase shifts, we further propose a suboptimal scheme leveraging the Newton’s method, which significantly reduces the computational burden while achieving performance close to the optimal SEE. Extensive numerical results validate the effectiveness of the proposed schemes, showing significant SEE improvements compared to traditional channel-capacity-based secure transmission scheme.

# Content of Code Package
This script performs the numerical evaluation of various RIS and beamforming schemes in a secure two-way relay (TWR) communication system. The results are accumulated over multiple iterations to analyze the overall performance. The following six schemes are implemented:

CCB (ccopt function): Implements the conventional channel capacity-based transmission scheme without employing physical layer security (PLS) techniques.Computes the secrecy rate (ccsr) and secrecy energy efficiency (cceta).

Random (ranopt function): A baseline scheme where both the RIS phase shifts and beamforming matrices are randomly selected from the feasible set. Computes the secrecy rate (ransr) and secrecy energy efficiency (raneta).

RandPhase (ranphiopt function): A scheme where the RIS phase shifts are randomly chosen, but the beamforming matrix is optimized iteratively. Computes the secrecy rate (ranphioptsr) and secrecy energy efficiency (ranphiopteta).

SR-lopt (ssropt function): A locally optimal secrecy rate (SR) maximization scheme. Jointly optimizes the RIS phase shifts and beamforming matrices to maximize the secrecy rate. Computes the secrecy rate (sroptsr), secrecy energy efficiency (sropteta), and optimized RIS phase shifts (sroptphi).

SEE-sub (seeoptsub function): A low-complexity suboptimal secrecy energy efficiency (SEE) maximization scheme. Computes the secrecy rate (seoptsubsr) and secrecy energy efficiency (seoptsubeta).

SEE-lopt (seeopt function): A locally optimal secrecy energy efficiency (SEE) maximization scheme. Jointly optimizes the RIS phase shifts and beamforming matrices for SEE enhancement. Computes the secrecy rate (seoptsr) and secrecy energy efficiency (seopteta).

For each scheme, the computed secrecy rate and secrecy energy efficiency values are accumulated to obtain the total performance metrics over multiple iterations.

# Requirement
To run this code, you need to install CVX: Matlab software for disciplined convex programming, version 2.1.

Download CVX: https://cvxr.com/cvx (Mar. 2014)

Ensure that CVX is properly installed and added to the MATLAB path before running the simulations. 
