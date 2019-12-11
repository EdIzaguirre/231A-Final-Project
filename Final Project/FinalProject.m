%% Rocket Landing with Grid Fins

% Constraints
Fmax = 1690*1000;
dmax = 5*pi/180;
umin = [0; -dmax*Fmax; -0.05; -0.05];
umax = [Fmax; dmax*Fmax; 0.05; 0.05];
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
u = sdpvar(4,N);
    
Q = diag([1 15 10 15])/(1e8);

constraints = [z(:,1) == z0, z(:,end) == zN];
cost = 0;

xBar = xBar();

for k = 1:N % stuff for all k to N-1
    xbar_k = xBar(k,:);
    cost = cost + z(:,k)'*Q*z(:,k) + u(2,k)^2 + (u(1,k)/Fmax)^2 + u(3,k)^2 + u(4,k)^2; %objective
    constraints = [constraints z(:,k+1) == RocketDynTrajectory(z(:,k),u(:,k),xbar_k) , umin<= u(:,k) <= umax]; %dynf and u constr
end

for k = 1:N+1 %stuff for all k to N
    constraints = [constraints zmin <= z(:,k)<= zmax]; %constr z_n+1
end

options = sdpsettings('solver','quadprog');
% options = sdpsettings('verbose',1,'solver','fmincon','usex0',1);
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
ylabel('\fin1 (degrees)')
plot(linspace(0,10,N), uOpt(4,:)./uOpt(1,:)*180/pi)
xlabel('t (s)')
ylabel('\fin2 (degrees)')











