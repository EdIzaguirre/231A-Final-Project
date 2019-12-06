clear all; clc;
load train_data_midterm
%% Load parameters
M = param.M;
g = param.g;
A = param.A;
B = param.B;
C = param.C;
%% Define state and input constraints
Fmax0 = param.mumax*M*g; umin = -Fmax0;
umax = Fmax0;
%% Grid the indipendent variable (same as gridding time) 
track_length = profile(end,1); % in meters
dp = 10; % sampling in space
p_sampled = 0:dp:track_length; %sampled train position
N_p = length(p_sampled);
%% Grid the state space
v_sampled = 0.1:0.1:max(maxspeed(p_sampled)); % grid from 0.1 to the max speed over the track a 
N_v = length(v_sampled);
v_idx_set = 1:N_v;
%% This funtion computes the next speed given the current speed, input and position
comp_v_next = @(v,u,p) v+dp/(v*M)*(-A-B*v-C*v^2-M*g*slope(p)-M*6/radius(p)+u);
%% This funtion computes the input to bring speed v to speed v next
% at position p
comp_u = @(v,v_next,p) M*(v_next-v)/(dp)*v-(-A-B*v-C*v^2-M*g*slope(p)-M*6/radius(p));
%% Define stage cost
Jstage = @(v,u) dp/v;
%% Initialization of Cost?to?Go 
for i = v_idx_set
if v_sampled(i) > maxspeed(p_sampled(N_p)) 
    J(N_p,i) = inf;
else
    J(N_p,i) = dp/v_sampled(i);
end
end
                             
 Jtogo = @(v) interpn(v_sampled,J(N_p,:),v,'linear');

 Jopt{N_p} = @(v) 0;
%% Perform Dynamic Programming


 tic
for p_indx = N_p-1:-1:1
fprintf('Solving DP at position = %i (meters) \n',p_sampled(p_indx)); 
J(p_indx,:) = inf(1,N_v);
u(p_indx,:) = nan(1,N_v);
p = p_sampled(p_indx);
for i = v_idx_set 
    v = v_sampled(i);
    Jbest = inf;
    Ubest = nan;
    um = comp_u(v,0.1,p);
    uM = comp_u(v,maxspeed(p+dp),p);
    for uGrid = linspace(max(umin,um),min(uM,umax),10)
        vnext = comp_v_next(v,uGrid,p);
        if vnext>=0.1 && vnext <= maxspeed(p+dp)
            Jactual = Jstage(v,uGrid) + Jopt{p_indx+1}(vnext); %Jtogo(vnext);
            if ~isnan(Jactual) 
                if Jactual< Jbest
                    Jbest = Jactual;
                    Ubest = uGrid;
                end
            end
        end
    end
    if v> maxspeed(p)
        Jbest = inf;
    Ubest = nan;
    end
    JoptArray(i,p_indx) = Jbest;
    UoptArray(i,p_indx) = Ubest;
end
Jopt{p_indx} = @(v) interpn(v_sampled,JoptArray(:,p_indx),v,'linear');
Uopt{p_indx} = @(v) interpn(v_sampled,UoptArray(:,p_indx),v,'linear');
end
fprintf('Total solution time: %i\n',toc);
%% Simulate the system from the initial condition v(0)=0.1 and p(0)=0 and plot

vOpt = zeros(1,N_p);
uOpt = zeros(1,N_p -1);
vOpt(1,1) = 0.1;%param.V0;

for i = 1:N_p-1
    uOpt(i) = Uopt{i}(vOpt(:,i));
    vOpt(:,i+1) = comp_v_next(vOpt(:,i),uOpt(:,i),p_sampled(i));
end

%%
subplot(3,1,1)
plot(p_sampled,maxspeed(p_sampled),p_sampled,vOpt)
legend('vmax','actual velocity')
subplot(3,1,2)
plot(p_sampled(1:end-1),uOpt)
xlabel('position')
ylabel('Optimal Control')
slopes = [];
subplot(3,1,3)
for p = p_sampled
    slopes = [slopes slope(p)];
end
plot(p_sampled,slopes)
xlabel('position')
ylabel('slope')


