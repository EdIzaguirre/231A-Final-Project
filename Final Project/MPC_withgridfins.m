%% Rocket Landing with Grid Fins (Linearized)
clc; clear;
% Constraints
%Fmax = 1690*1000;
Fmax = 856000;
dmax = 5*pi/180;
umin = [0; -dmax*Fmax; 0; 0];
umax = [Fmax; dmax*Fmax; 90*pi/180; 90*pi/180];  % 200 is miniumum for u3 and u4
tmin = -20*pi/180;
tmax = 20*pi/180;
zmin = [tmin; -100; -3000; -100]; zmax = [tmax; 100; 0; 500];

% Initial conditions
v0 = 110;
alt0 = -1410;
t0 = 10*pi/180;
z0 = [t0; 0; alt0; v0];
zN = [0;0;0;0];

% % LQR
% A = [1 0.1 0 0; 0 1 0 0; 0 0 1 0.1; 0 0 0 1]; %Abar and Bbar from z at N+1
% B = [0 0 0 0; 0 -0.3465 0 0; 0 0 0 0; 0 0 -0.0004 -0.0004]; 
% Q = diag([1 15 10 15])/(1e8);
Q = diag([0.1 0.1 0.1 0.1]);
% R = 1;
% [Finf,P]=dlqr(A,B,Q,R);
% AclF=A-B*Finf;
% nu = size(B,2);
% 
% % constraint sets represented as polyhedra
% Z = Polyhedron('lb',zmin,'ub',zmax);
% U = Polyhedron('lb.',umin,'ub',umax);
% % remember to convert input constraints in state constraints
% S=Z.intersect(Polyhedron('H',[-U.H(:,1:nu)*Finf U.H(:,nu+1)]));
% tic
% OinfF=max_pos_inv(AclF,S);
% toc
% Af=OinfF.H(:,1:end-1);
% bf=OinfF.H(:,end);

% Define sampling time
TS = 0.2525;


% Define horizon
N = 100;
M = 130;
zOpt = zeros(4,M+1);
zOpt(:,1) = z0;
uOpt = zeros(4,M);

xBar = xBar();

alpha1 = 1/((Fmax)^2);
alpha2 = 0.001;
alpha3 = 0.001;
alpha4 = 0.001;

% Creating noise with mean 0 and variance 10
mu=0; sigma=sqrt(10);
rng(0);
noise = sigma*randn(4,N)+mu;

%% Optimization
tic
for t = 1:M
%     xBar = xBar1(t:t+N,:);
    disp(t);
    z = sdpvar(4,N+1);
    u = sdpvar(4,N);

    constraints = [z(:,1) == zOpt(:,t) z(:,end) == zN];
%     constraints = [constraints, Af*z(:,N+1)<=bf];
    cost = 0;
%     cost = cost + z(:,N+1)'*P*z(:,N+1);

    for k = 1:N % stuff for all k to N-1
        xbar_k = xBar(k,:);
        xbar_kNext = xBar(k+1,:);
        cost = cost + z(:,k)'*Q*z(:,k) + alpha1*u(1,k)^2 + alpha2*u(2,k)^2 + alpha3*u(3,k)^2 + alpha4*u(4,k)^2;
%         cost = cost + z(:,k)'*Q*z(:,k) + (u(1,k)/Fmax)^2 + u(2,k)^2 + u(3,k)^2 + u(4,k)^2; %objective
        constraints = [constraints z(:,k+1) == RocketDynTrajectory(z(:,k),u(:,k),xbar_k,xbar_kNext, noise(:,k)) , umin<= u(:,k) <= umax];
        constraints = [constraints zmin <= z(:,k)<= zmax];
    end

    if k == N %stuff for all k to N
        constraints = [constraints zmin <= z(:,N+1)<= zmax]; %constr z_n+1
    end

    options = sdpsettings('solver','quadprog');
    % options = sdpsettings('verbose',1,'solver','fmincon','usex0',1);
    sol = optimize(constraints, cost, options);
    if sol.problem == 0
        zO = value(z);
        uO = value(u);
    else
        feas = 0;
        disp("infeasible solution");
        disp("exiting the simulation");
%         zO = [];
%         uO = [];
        return;
    end
    
    uOpt(:,t) = uO(:,1);
    if t>100
        k =100;
    else
        k = t;
    end
    zOpt(:,t+1) = RocketDynTrajectory(zOpt(:,t),uOpt(:,t),xbar_k,xbar_kNext, noise(:,k));
end
toc
%% Plotting
figure;
subplot(2,2,1)
plot(linspace(0,10,M+1), zOpt(1,:))
xlabel('t')
ylabel('\theta (rad)')
subplot(2,2,2)
plot(linspace(0,10,M+1), zOpt(2,:))
xlabel('t')
ylabel('\omega (rad/s)')
subplot(2,2,3)
plot(linspace(0,10,M+1), zOpt(3,:))
xlabel('t (s)')
ylabel('h (m)')
subplot(2,2,4)
plot(linspace(0,10,M+1), zOpt(4,:))
xlabel('t (s)')
ylabel('v (m/s)')

figure;
subplot(1,4,1)
plot(linspace(0,10,M), uOpt(1,:))
xlabel('t (s)')
ylabel('F (N)')
subplot(1,4,2)
plot(linspace(0,10,M), uOpt(2,:)*180/pi)
xlabel('t (s)')
ylabel('\delta (degrees)')
subplot(1,4,3)
plot(linspace(0,10,M), uOpt(3,:)*180/pi)
xlabel('t (s)')
ylabel('fin1 (degrees)')
subplot(1,4,4)
plot(linspace(0,10,M), uOpt(4,:)*180/pi)
xlabel('t (s)')
ylabel('fin2 (degrees)')





