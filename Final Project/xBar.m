function [xBar] = xBar()
    %% Finding fixed trajectory
    % Note: Borelli said to set all grid fin inputs to 0 (u3 and u4), set the
    % gimbal angle to 0 (u2) and set F (u1) to some function slightly greater
    % than gravity, in order to bring the system to rest. I will simulate as a
    % constant force/constant acceleration problem.

    %% Setting up super simple rocket problem

    % Parameters
    m = 27648;
    l = 70;
    J = 1/16*m*l^2;
    g = 9.8;
    
    % Initial conditions
    deltaT = 10; %seconds
    deltaH = 3000; % meters
    v0 = 205.2; % meters

    % Finding constant force to get rocket to land
    F = m * (g - 2*(deltaH - v0*deltaT)/(deltaT^2));

    %% Simulating for N steps

    % Define sampling time
    TS = 0.1; h0 = -3000;

    % Define horizon
    N = deltaT/TS;  % 100 Samples

    thetaList = zeros(N,1); omegaList = zeros(N,1); hListGood = h0; vList = v0;

    for i = 1:N-1
        nextHeight = h0 + 1/2*(g - F/m)*(TS)^2 + v0*TS;
        nextVelocity = (g - F/m)*TS + v0;

        % Append to lists and reinitializ v0 and h
        hListGood = [hListGood; nextHeight]; vList = [vList; nextVelocity]; 
        v0 = vList(end); h0 = hListGood(end);
    end

    %% Finding fixed trajectory
    % Note: Borelli said to set all grid fin inputs to 0 (u3 and u4), set the
    % gimbal angle to 0 (u2) and set F (u1) to some function slightly greater
    % than gravity, in order to bring the system to rest. I will simulate as a
    % constant force/constant acceleration problem.

    %% Setting up super simple rocket problem

    % Initial conditions
    deltaT = 10; %seconds
    v0 = -205.2; % meters

    % Finding constant force to get rocket to land
    F = m*g - m*v0/deltaT

    %% Simulating for N steps

    % Define sampling time
    TS = 0.1; h0 = 3000;

    % Define horizon
    N = deltaT/TS;  % 100 Samples

    thetaList = zeros(N,1); omegaList = zeros(N,1); hList = h0; vListGood = v0;

    for i = 1:N-1
        nextVelocity = (F/m - g)*TS + v0;
        % Append to lists and reinitializ v0
        vListGood = [vListGood; nextVelocity]; 
        v0 = vListGood(end);
    end

    vListGood = -vListGood;

    % Compile
    
    
    thetaList = zeros(N,1); omegaList = zeros(N,1);
    hListGood; vListGood;
    u1 = F*ones(N,1);
    u2 = zeros(N,1);
    u3 = zeros(N,1);
    u4 = zeros(N,1);
    
    xBar = [thetaList, omegaList, hListGood, vListGood, u1, u2, u3, u4];
end





















