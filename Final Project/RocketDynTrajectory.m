function znext = RocketDynTrajectory(z,u,xbar,xbarNext)
TS = 0.1;
m = 27648;
l = 70;
J = 1/16*m*l^2;
g = 9.8;
rho = 1.225;
Cd = 0.14;
A_tot = 1.5;
b = 9;
Fd1 = 1/2*rho*z(4)^2*Cd*A_tot*sin(u(3));
Fd2 = 1/2*rho*z(4)^2*Cd*A_tot*sin(u(4));

syms theta omega h v u1 u2 u3 u4
Syms = [theta omega h v u1 u2 u3 u4];
    
% Getting Linear Dynamics and Placing into Functions
[Ac,Bc] = linDynamics(); 
Ad = eye(4) + TS*Ac;
Bd = TS*Bc;
A_bar = double(subs(Ad,Syms(1:8),xbar(1:8)));
B_bar = double(subs(Bd, Syms(1:8),xbar(1:8)));

deltaX = z-xbar(1:4)';
deltaU = u-xbar(5:8)';
deltaXnext = A_bar*deltaX + B_bar*deltaU;
znext = xbarNext(1:4)' + deltaXnext;

end