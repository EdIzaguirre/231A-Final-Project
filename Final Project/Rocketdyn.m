function znext = Rocketdyn(z,u)
TS = 0.1;
m = 27648;
l = 70;
J = 1/16*m*l^2;
g = 9.8;
rho = 1.225;
Cd = 0.14;
A_tot = 1.5;

znext(1,1) = z(1) + TS*z(2);
znext(2,1) = z(2) - TS*l/2/J*u(2);
znext(3,1) = z(3) + TS*z(4);
znext(4,1) = z(4) +TS*(g - l/m*u(1) - 4*1/2*rho*z(4)^2*Cd*A_tot*sin(u3)/m);
end