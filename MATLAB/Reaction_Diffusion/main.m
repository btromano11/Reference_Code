%Solves the below Reaction-Diffusion system

%Ut = lambda(A)*U - omega(A)*V + D1*(Lap)^2*U
%Vt = omega(A)*U + lambda(A)*V + D2*(Lap)^2*V

%where:
    %A^2 = U^2 + V^2
    %(Lap)^2 = dx^2 + dy^2
    %lambda(A) = 1 - A^2
    %omega(A) = 1 - beta*A^2

clear all;close all;clc

solvec = spec_solve();

solvec1 = real(solvec);
solvec2 = imag(solvec);
save A1.dat solvec1 -ascii
save A2.dat solvec2 -ascii

solvec = cheb_solve();
save A3.dat solvec -ascii
