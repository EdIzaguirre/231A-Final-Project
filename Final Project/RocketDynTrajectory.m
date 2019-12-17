function znext = RocketDynTrajectory(z,u,xbar,xbarNext,noise)
    syms theta omega h v u1 u2 u3 u4
    Syms = [theta omega h v u1 u2 u3 u4];


    D = [0.001;
        0.001;
        0.1;
        0.1];   % Angle noise needs to be smaller (due to radians) than other noise.
    
    TS = 0.1;

    % Getting Linear Dynamics and Placing into Functions
    [Ac,Bc] = linDynamics(); 
    Ad = eye(4) + TS*Ac;
    Bd = TS*Bc;
    A_bar = double(subs(Ad,Syms(1:8),xbar(1:8)));
    B_bar = double(subs(Bd, Syms(1:8),xbar(1:8)));

    deltaX = z-xbar(1:4)';
    deltaU = u-xbar(5:8)';  
    deltaXnext = A_bar*deltaX + B_bar*deltaU + D.*noise;
    znext = xbarNext(1:4)' + deltaXnext;

end