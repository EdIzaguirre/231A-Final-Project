function A = linDynamics()
    clear; clc;

    %% Begin with continous time dynamics

    syms theta omega h v u1 u2 u3 u4

    Syms = [theta omega h v u1 u2 u3 u4];
    % Parameters
    m = 27648;
    l = 70;
    J = 1/16*m*l^2;
    g = 9.8;

    % Drag constants/variables
    rho = 1.225;
    Cd = 0.14;
    A_tot = 1.5;
    b = 9;
    Fd1 = 1/2*rho*v^2*Cd*A_tot*sin(u3);
    Fd2 = 1/2*rho*v^2*Cd*A_tot*sin(u4);

    % Definition of functions
    f1 = omega;
    f2 = -l/(2*J)*u1*u2 + Fd1*b/(2*J) - Fd2*b/(2*J);
    f3 = v;
    f4 = g - (1/m)*u1 - (1/m)*(Fd1 + Fd2);

    % Obtain linearized dynamics 
    A = jacobian([f1, f2, f3, f4],Syms);

end