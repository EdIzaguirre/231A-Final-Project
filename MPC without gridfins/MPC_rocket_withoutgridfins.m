clear
% Parameters
m = 27648;
l = 70;
J = 1/16*m*l^2;
g = 9.8;

% Constraints
Fmax = 1690*1000;
dmax = 5*pi/180;
umin = [0; -dmax*Fmax];
umax = [Fmax; dmax*Fmax];
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
N = 55;
M = 100;
zOpt = zeros(4,M+1);
zOpt(:,1) = z0;
uOpt = zeros(2,M);


%%
for t = 1:M
    disp(t);
    z = sdpvar(4,N+1);
    u = sdpvar(2,N);

    Q = diag([1 1 5 5])/(1e8);

    constraints = [z(:,1) == z0, z(:,end) == zN];
    cost = 0;
    for k = 1:N
        cost = cost + z(:,k)'*Q*z(:,k) + u(2,k)^2 + (u(1,k)/Fmax)^2;
        constraints = [constraints z(:,k+1) == Rocketdyn(z(:,k),u(:,k)), umin<= u(:,k) <= umax];
    end

    for k = 1:N+1
        constraints = [constraints zmin <= z(:,k)<= zmax];
    end

    options = sdpsettings('verbose',1,'solver','quadprog');
    sol = optimize(constraints, cost, options);
    
    if sol.problem == 0
        zOpt = value(z);
        uOpt = value(u);
    else
        feas = 0;
        disp("infeasible solution");
        disp("exiting the simulation");
        xOpt = [];
        uOpt = [];
        return;
    end
    
    uOpt(:,t) = u(:,1);
    zOpt(:,t+1) = Rocketdyn(zOpt(:,t),uOpt(:,t));
end
    

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