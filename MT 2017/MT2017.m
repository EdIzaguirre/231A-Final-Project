%% Problem 1(b)
x1 = sdpvar(1,1);x2 = sdpvar(1,1);
cost = -(x1^2+0.5*x2^2)/2;
constraints = [x1+2*x2 == 1, 0<= x1<=1, 0<=x2<=1];
options = sdpsettings('verbose', 0, 'solver', 'quadprog');
optimize(constraints, cost, options);
x1Opt = double(x1)
x2Opt = double(x2)
JOpt = double(cost)
v = dual(constraints(1))
u1 = dual(constraints(2))
u2 = dual(constraints(3))

%% Derek Practice 1b
A = [-1 0 ; 1 0 ; 0 -1; 0 1];
b = [0; 1; 0; 1];
cost = @(x) -(x1^2+0.5*x2^2)/2;
[x,fval,exitflag,output,lambda] = fmincon(cost,A,b);
u = lambda.ineqlin

%% Problem 1(b) (ii)
g_f = [-x1Opt;-x2Opt/2];
g_h = [1;2]; % Derivatives of dyanmic constraint
g_g1 = [-1;0];
g_g2 = [1;0];
g_g3 = [0;-1];
g_g4 = [0;1];


epsilon = 1e-8;
check1 = all(abs(g_f + u1(1)*g_g1 + u1(2)*g_g2 + u2(1)*g_g3 + u2(2)*g_g4 + v*g_h)<epsilon);
check2 = abs(u1(1)*x1Opt)<epsilon && abs(u1(2)*(x1Opt-1))<epsilon && abs(u2(1)*x2Opt)<epsilon && abs(u2(2)*(x2Opt-1))<epsilon;
check3 = u1(1) >=0 && u1(2) >=0 && u2(1) >=0 && u2(2) >=0;
check4 = abs(x1Opt+2*x2Opt-1)<epsilon; % h_j >=0;
check5 = x1Opt>= 0 && x1Opt<= 1 && x2Opt>= 0 && x2Opt<= 1; %g_i(zstar) >=0
KKTsat = check1 && check2 && check3 && check4 && check5;

%% Problem 1(c) ii
c = [0;1];
A = [3 -1;
     -3 -1;
     -1 0;
     1 0];
b = [-2;2;3;10];
[x,fval] = linprog(c,A,b)
z = x(1);

%% Problem 3
clear
% Parameters
    m = 27648;
    l = 70;
    J = 1/16*m*l^2;
    g = 9.8;
    % Constraints
Fmax = 1690*1000;
dmax = 5*pi/180;
umin = [0; -dmax*Fmax; 0];
umax = [Fmax; dmax*Fmax; .17*Fmax];
tmin = -20*pi/180;
tmax = 20*pi/180;
zmin = [tmin; -100; -3000; -100]; zmax = [tmax; 100; 0; 500];
    % Initial conditions
v0 = 205.2;
alt0 = -1228; %negative because the Z?axis is positive pointing downward 
t0 = 10*pi/180;
z0 = [t0; 0; alt0; v0];
zN = [0;0;0;0];
    % Define sampling time
TS = 0.1;
    % Define horizon
N = 10/TS;

%%
% ok to seperate the terms because the cost is all within the same k
tic
z = sdpvar(4,N+1);
u = sdpvar(3,N);

Q = diag([1 1 5 5])/(1e8);

constraints = [z(:,1) == z0, z(:,end) == zN];
cost = 0;

for k = 1:N % stuff for all k to N-1
    cost = cost + z(:,k)'*Q*z(:,k) + u(2,k)^2 + (u(1,k)/Fmax)^2 + u(3,k)^2; %objective
    constraints = [constraints z(:,k+1) == Rocketdyn(z(:,k),u(:,k)) , umin<= u(:,k) <= umax]; %dynf and u constr
end

for k = 1:N+1 %stuff for all k to N
    constraints = [constraints zmin <= z(:,k)<= zmax]; %constr z_n+1
end

options = sdpsettings('verbose',1,'solver','fmincon','usex0',1);
diagnostics = optimize(constraints, cost, options);
zOpt = value(z)
uOpt = value(u)
toc
%%
figure;
subplot(2,2,1)
plot(linspace(0,10,N+1), zOpt(1,:))
xlabel('t')
ylabel('\theta (rad)')
subplot(2,2,2)
plot(linspace(0,10,N+1), zOpt(2,:))
xlabel('t')
ylabel('\omega (rad/s)')
subplot(2,2,3)
plot(linspace(0,10,N+1), zOpt(3,:))
xlabel('t (s)')
ylabel('h (m)')
subplot(2,2,4)
plot(linspace(0,10,N+1), zOpt(4,:))
xlabel('t (s)')
ylabel('v (m/s)')

figure;
subplot(1,3,1)
plot(linspace(0,10,N), uOpt(1,:))
xlabel('t (s)')
ylabel('F (N)')
subplot(1,3,2)
plot(linspace(0,10,N), uOpt(2,:)./uOpt(1,:)*180/pi)
xlabel('t (s)')
ylabel('\delta (degrees)')
subplot(1,3,3)
plot(linspace(0,10,N), uOpt(3,:)./uOpt(1,:)*180/pi)
xlabel('t (s)')
ylabel('\delta (degrees)')

%% Problem 3(c)
N = 5/TS;
z = sdpvar(4,N+1);
u = sdpvar(2,N);

Q = diag([1 1 5 5]);

constraints = [z(:,1) == z0, z(:,end) == zN];
cost = 0;
for k = 1:N
    cost = cost + z(:,k)'*Q*z(:,k) + u(2,k)^2 + (u(1,k)/Fmax)^2;
    constraints = [constraints z(:,k+1) == Rocketdyn(z(:,k),u(:,k)), umin<= u(:,k) <= umax];
end

for k = 1:N+1
    constraints = [constraints zmin <= z(:,k)<= zmax];
end

options = sdpsettings('verbose',1,'solver','fmincon','usex0',z0);
diagnostics = optimize(constraints, cost, options);
zOpt = value(z);
uOpt = value(u);

%%
figure;
subplot(2,2,1)
plot(linspace(0,10,N+1), zOpt(1,:))
xlabel('t')
ylabel('\theta (rad)')
subplot(2,2,2)
plot(linspace(0,10,N+1), zOpt(2,:))
xlabel('t')
ylabel('\omega (rad/s)')
subplot(2,2,3)
plot(linspace(0,10,N+1), zOpt(3,:))
xlabel('t (s)')
ylabel('h (m)')
subplot(2,2,4)
plot(linspace(0,10,N+1), zOpt(4,:))
xlabel('t (s)')
ylabel('v (m/s)')

figure;
subplot(1,2,1)
plot(linspace(0,10,N), uOpt(1,:))
xlabel('t (s)')
ylabel('F (N)')
subplot(1,2,2)
plot(linspace(0,10,N), uOpt(2,:)./uOpt(1,:)*180/pi)
xlabel('t (s)')
ylabel('\delta (degrees)')





























