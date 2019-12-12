function znext = RocketDynTrajectory(z,u,xbar)
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
A = linDynamics(); 
A_bar = double(subs(A,Syms,xbar));
new_dyn = A_bar*transpose(Syms);
f_1 = new_dyn(1); f_2 = new_dyn(2); f_3 = new_dyn(3); f_4 = new_dyn(4);

%% Converting from Sym to YALMIP
s = char(f_1);
s = strrep(s,char(omega),'z(2)');
f_1 = eval(s);

s = char(f_2);
s = strrep(s,char(u1),'u(1)');
s = strrep(s,char(u2),'u(2)');
s = strrep(s,char(u3),'u(3)');
s = strrep(s,char(u4),'u(4)');
f_2 = eval(s);

s = char(f_3);
s = strrep(s,char(v),'z(4)');
f_3 = eval(s);

s = char(f_4);
s = strrep(s,char(u1),'u(1)');
s = strrep(s,char(u3),'u(3)');
s = strrep(s,char(u4),'u(4)');
f_4 = eval(s);

%% Propogating Dynamics
znext(1,1) = z(1) + TS*f_1;
znext(2,1) = z(2) + TS*f_2;
znext(3,1) = z(3) + TS*f_3;
znext(4,1) = z(4) + TS*f_4;

end