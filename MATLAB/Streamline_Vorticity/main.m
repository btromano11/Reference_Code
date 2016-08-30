%Solves the below Streamline-Vorticity System

%wt + [phi,w] = v(Lap)^2w

%where:
    %[phi,w] = phi_x*w_x - phy_y*w_y
    %(Lap)^2 = dx^2 + dy^2
    %Lap)^2phi = w


clear variables;close all;clc;

L=10;
n = 64;

%%%%%FFT
y = solvefft(L,n);
save A1.dat y -ascii
waterfall(y);







