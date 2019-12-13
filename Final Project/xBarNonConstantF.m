function [xBar] = xBarNonConstantF()
    clc; clear;
    %% Finding new function for F
%     syms a b
% 
%     func1 = 1690000 == (1/a + b);
%     func2 = 0 == 1/(10+a) + b;
%     [asol, bsol] = solve(func1,func2);
% 
%     double(asol);
%     double(bsol);
% 
%     % Didn't work
%     %% Checking momentum loss
% 
%     % The equation f(t) = -169000*t + 1690000 seems to be able to roughly get
%     % what we want out of our force (in 10 seconds).
% 
%     % Checking desired change in momentum
%     m = 27648;
%     l = 70;
%     J = 1/16*m*l^2;
%     g = 9.8;
%     v0 = 205.2;
%     alt0 = -1228;   % OG height
% 
%     initMomentum = m*v0;
%     finalMomentum = 0;
% 
%     desiredMomentumLoss = finalMomentum - initMomentum % -5.6734*10^6 N*s of impulse needed by rocket thrust.
% 
%     % Calculating actual change in momentum
%     tmin = 0;
%     tmax = 10;
%     fun = @(t) -169000*t + 1690000;
% 
%     actualMomentumLoss = integral(fun,tmin,tmax)
% 
%     abs(actualMomentumLoss) - abs(desiredMomentumLoss)
% 
%     % Not enough momentum taken out

    %% Via Batch Optimization

    % Constraints
    Fmax = 1690*1000;  
    dmax = 5*pi/180;    % 5 Degrees
    umin = [0; 0; 0; 0];
    umax = [Fmax; 0; 0; 0]; 
    tmin = 0;
    tmax = 0;   %20 degrees
    zmin = [0; 0; -3000; -100]; zmax = [0; 0; 0; 500];

    % Initial conditions
    v0 = 205.2;
    alt0 = -1228;
    t0 = 0; % 0 degrees
    z0 = [t0; 0; alt0; v0];
    zN = [0;0;0;0];
    u0 = [Fmax;0;0;0];
    uN = [0;0;0;0];

    % Define sampling time
    TS = 0.1;

    % Define horizon
    N = 10/TS;

    %% Optimization
    tic
    z = sdpvar(4,N+1);
    u = sdpvar(4,N);

    Q = diag([1 15 10 15])/(1e8);

    constraints = [z(:,1) == z0, z(:,end) == zN u(:,1) == u0];
    cost = 0;

    for k = 1:N % stuff for all k to N-1
        if k == 1
            cost = cost + z(:,k)'*Q*z(:,k) + (u(1,k)/Fmax)^2 + u(2,k)^2 + u(3,k)^2 + u(4,k)^2; %objective
        else
            cost = cost + z(:,k)'*Q*z(:,k) + (u(1,k)/Fmax)^2 + u(2,k)^2 + u(3,k)^2 + u(4,k)^2 + norm(u(1,k)-u(1,k-1))^2; %objective
        end
        constraints = [constraints z(:,k+1) == simplifiedRocketDyn(z(:,k),u(:,k)) , umin<= u(:,k) <= umax]; %dynf and u constr
    end

    for k = 1:N+1 %stuff for all k to N
        constraints = [constraints zmin <= z(:,k)<= zmax]; %constr z_n+1
    end

    options = sdpsettings('solver','quadprog');
    % options = sdpsettings('verbose',1,'solver','fmincon','usex0',1);
    diagnostics = optimize(constraints, cost, options);
    zOpt = value(z);
    uOpt = value(u);
    toc
    
    % Compile all vectors for xBar
    thetaList = zeros(N,1); 
    omegaList = zeros(N,1);
    hList = zOpt(3,1:N)';
    vList = zOpt(4,1:N)';
    u1List = uOpt(1,:)';
    u2List = zeros(N,1);
    u3List = zeros(N,1);
    u4List = zeros(N,1);
    xBar = [thetaList, omegaList, hList, vList, u1List, u2List, u3List, u4List];
    
    %% Plotting
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
    subplot(1,4,1)
    plot(linspace(0,10,N), uOpt(1,:))
    xlabel('t (s)')
    ylabel('F (N)')
    subplot(1,4,2)
    plot(linspace(0,10,N), uOpt(2,:)*180/pi)
    xlabel('t (s)')
    ylabel('\delta (degrees)')
    subplot(1,4,3)
    plot(linspace(0,10,N), uOpt(3,:)*180/pi)
    xlabel('t (s)')
    ylabel('fin1 (degrees)')
    subplot(1,4,4)
    plot(linspace(0,10,N), uOpt(4,:)*180/pi)
    xlabel('t (s)')
    ylabel('fin2 (degrees)')

end
