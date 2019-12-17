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
    g = 9.8;

    % Initial conditions
    deltaT = 10; %seconds
    v0 = -205.2; % meters/sec
    y0 = 1061; % meters

    % Finding constant force to get rocket to land
    F = m * (g + (0 - (v0))/deltaT);

    %% Simulating for N steps

    % Calculate Acceleration
    a = F/m - g; %m/s^2

    % Getting Proper Time Spacing 
    N = 101;  % 100 Samples
    
    totTime = linspace(0,10,N);
    TS = totTime(2) - totTime(1);  

    % Begin Simulation
    vList = v0; hList = y0; t = 0;

    for i = 2:N
        nextVelocity = a*TS + v0;
        nextHeight = y0 + (1/2) * a * TS^2 + v0 * TS;
        t = t+TS;
        
        % Append to lists and reinitialize v0 and h
        vList = [vList; nextVelocity]; hList = [hList; nextHeight]; 
        v0 = vList(end); y0 = hList(end);
    end

    % Compile all vectors for xBar
    thetaList = zeros(N,1); omegaList = zeros(N,1);
    u1List = F*ones(N,1);
    u2List = zeros(N,1);
    u3List = zeros(N,1);
    u4List = zeros(N,1);
    xBar = [thetaList, omegaList, -hList, -vList, u1List, u2List, u3List, u4List];
    % Negative because I solved the above problems assuming up is positive.
    
    %% Checking to see if the trajectory looks good
%     figure();
%     subplot(2,1,1)
%     plot(linspace(1,10,101),-hList)
%     title('Height vs. Time')
%     xlabel('Time (sec)')
%     ylabel('Height (m) (Negative is up)')
%     hold on;
%     subplot(2,1,2)
%     plot(linspace(1,10,101),-vList)
%     title('Velocity vs. Time')
%     xlabel('Time (sec)')
%     ylabel('Velocity (m/s) (Positive is down)')

end





















