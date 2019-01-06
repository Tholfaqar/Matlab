%% 220B Final Project
%% Charging Controller

yalmip('clear')
clear all
close all
%% Model Declaration and Parameter Choice

% In this section, we present a 9-section main feeder line model that is used in Example 1 of the paper.
% To see this model, refer to http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4111450


% | -------| -----| ..... -----| ------- |
%b1       b2      b3           b9       substation
%% Declaration of Bus Data (Real and Reactive Power Consumption at each Bus)

%first we have constant power draw. We will revise it to be a 24-hr periodic profile afterwards

%real power in kW, 1640kW is the real power consumption at bus 1 for example
RealPowerLoad = [1640 980 1150 780 1610 1598 1790 980 1840]';

%reactive power in kvar
ReactivePowerLoad = [200 130 60 110 600 1840 446 340 460]';
BusData = [RealPowerLoad,ReactivePowerLoad];
numBus = length(RealPowerLoad);

%% Declaration of Line Data

% lineMatrix(i,j) contains the impedance from bus i to bus j. It will be in the form r+ix where r is the resistance
% and x is the reactance.
% bus 10 is the substation and the voltage there is 23kV
lineMatrix = zeros(10,10);
lineMatrix(1,2) = 5.3434 + 1i*3.0264;
lineMatrix(2,3) = 4.7953 + 1i*2.7160;
lineMatrix(3,4) = 2.0552 + 1i*1.1640;
lineMatrix(4,5) = 0.9053 + 1i*0.7886;
lineMatrix(5,6) = 1.9831 + 1i*1.7276;
lineMatrix(6,7) = 0.6984 + 1i*0.6084;
lineMatrix(7,8) = 0.7463 + 1i*1.2050;
lineMatrix(8,9) = 0.0140 + 1i*0.6051;
lineMatrix(9,10) = 0.1233 + 1i*0.4127;

resistanceMatrix = real(lineMatrix);
reactanceMatrix = imag(lineMatrix);
%% Construction of R,X matrices

% We wish the express the bus voltage as a function of the substation voltage (assume constant), impedances of the line 
% and the load profiles at each bus. The general expression is given by:
% Vi = V0 - R(pl + pv) - X(ql - qg)

Vsub = 23e3;
%R(i,j) = sum of Rs common to Bus i and Bus j. For example, R(1,2) = Sum of Rs from Bus 2 to substation.
R = zeros(9,9);
tempR = 0;
for i = 0:8
    tempR = tempR + resistanceMatrix(9-i,10-i);
    R(1:end-i,9-i) = tempR * ones(9-i,1);
end

X = zeros(9,9);
tempX = 0;
for j = 0:8
    tempX = tempX + reactanceMatrix(9-j,10-j);
    X(1:end-j,9-j) = tempX*ones(9-j,1);
end

%% V_tilda

%V_tilda = Vsub - R*pload - X*qload (grouping all known terms together)
V_tilda = Vsub/Vsub - R*RealPowerLoad*10^3/(Vsub^2) - X*ReactivePowerLoad*10^3/(Vsub^2);
    
%V_tilda is in per unit. 

%% define constraints

% SOC constraits
% SOC_min <= SOC_des - SOC_cur <= SOC_max

% Reactive generation constraints
% qmin <= qg <= qmax

% Charging power constraint
% 0 <= u <= cmax

% Line voltage constraint
% Vmin - V_tilda <= y - Vnom <= Vmax - V_tilda


%set cnom to be 1/10 of max(RPL). Increasing this allows faster charging. If the design tracjectory is always cmax,
%then we should increase this value to make the problem more interesting. 
cnom = 10*max(RealPowerLoad)*1e3; 
cmax = 0.6*cnom;

% assuming a full charge (from 0% to 100%) takes two hours to complete. 
% SOC_nom = 120 mins * cnom. Increasing this makes charging 'slower' to complete. 
% We can reduce this to make simulation shorter.
SOC_nom = 120 * 60 * cmax;
SOC_min = 0.2 * SOC_nom;
SOC_max = 1.0 * SOC_nom;

SOC_maxpu = 1;

%This depends on physical constraint eg how many inverters connected and their capabilities. 
%Decreasing this makes the problem harder. (having less reactive power injection capability)
qnom = 2*1000*1e3; 
qmin = -1.0*qnom;
qmax = 1.0*qnom;

%all voltages are scaled to per unit
Vmin = -0.2;
Vmax = 0.2;
Vnom = 1;
%% state space definitions

% x (states) are the SOC_des - SOC_cur  of the EVs connected to the grid
% control inputs are the n reactive power generation as well as the charging rate 
% A = I, B = [0 -T], C = 0, D = [X -RK] where K is the EV connection matrix 
% In this section, I will put forward the SS realization as well as the constraints

% Definition of K, the EV connection matrix. K is a n*m matrix where n is the number of buses in the network
% and m is the number of EVs present in the network. K(i,j) == 1 means EV j is connected to bus i. In this example, 
% we assume constant number of EVs connected to the network. 
% Each EV can be connected to only one bus hence the column sum of the matrix is 1. 


% assume we have 9 EVs in the network

numEVs = 9;

% Different Scenarios of EV charging
% 1. one EV connected to each bus
K = eye(9);
% 2. all EVs connected to bus 1 (terminal bus)
% K = zeros(9,9);
% K(1,:) = ones(1,9); 
% 3. all EVs connected to bus 9 (nearest bus but with the heaviest load)
%K = zeros(9,9);
%K(9,:) = ones(1,9);
% 4. all EVs connected to the bus with the lightest load (bus 4)
%K = zeros(9,9);
%K(4,:) = ones(1,9);


% define SOC initial conditions

SOC_des = [1.0 0.8 0.9 0.8 1.0 0.8 0.6 0.7 1.0]';
SOC_start = [0.5 0.3 0.5 0.2 0.1 0.2 0.1 0.1 0.4]';
SOC_initial_conditions = SOC_des - SOC_start;

T = 5; %sampling time is 5 mins
A = eye(9);
B = [zeros(1,numBus) -T*60*ones(1,numEVs)/SOC_nom];
C = 0;
D = [X -R*K]/(Vsub^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation and Model Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nx = length(A); % Number of states
nuq = size(D,1); % Number of inputs
nu2 = nx; % Number of inputs
ny = size(D,1);  % number of outputs

%N = 24; %2 hours with 5 min steps
N = 4;

Nr = 24*2; % if 5 min intervals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stage 1: Periodic Load 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Construct a load profile that is 24 hour periodic. Given that our sampling time is 5min, the load profile
%is 24*60/5 = 288 long for each bus instead of a constant value

%Let's have a isosceles triangle profile (peaking in the middle of the day).Base multiplier is 0.7. Peak
%multiplier is 1.3.
a = 0.65;
b = 2.0;
r = a + (b-a).*rand(1,Nr/2);
r1 = sort(r,'ascend');
r2 = sort(r,'descend');
ref = zeros(1,Nr);
ref(1:Nr/2) = r1;
ref(Nr/2+1:end) = r2;
profile = ref;

RealPowerProfile = RealPowerLoad*profile;
ReactivePowerProfile = ReactivePowerLoad*profile;

figure
plot(profile)
title('Periodic Load Profile (p.u)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stage 1 V_tildaProfile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_tildaProfile = Vsub/Vsub - R*RealPowerProfile*10^3/(Vsub^2) - X*ReactivePowerProfile*10^3/(Vsub^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stage 1 MPC Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T1 = eye(ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stage 1 Define Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uqtil = sdpvar(repmat(nuq,1,Nr),repmat(1,1,Nr));
ytil = sdpvar(repmat(ny,1,Nr),repmat(1,1,Nr));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stage 1 Initialize  Control Desing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xref = [];
u2ref = [];
uqref = [];
yref = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stage 1 Setup Optimzation Problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

constraints =[];
for i = 1:Nr   
    constraints = [constraints, ytil{i} == D(:,1:nuq)*uqtil{i}];
    constraints = [constraints, qmin<= uqtil{i} <= qmax];   
    constraints = [constraints, Vmin*ones(ny,1)-V_tildaProfile(:,i) <= ytil{i}-Vnom*ones(ny,1) <= Vmax*ones(ny,1)-V_tildaProfile(:,i)]; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stage 1 Solve MPC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
objective = 0;
r = Vnom*ones(size(V_tildaProfile)) - V_tildaProfile;
for i = 1:Nr
    objective = objective + (ytil{i}-r(:,i))'*T1*(ytil{i}-r(:,i)); 
                
end
     

diag=optimize(constraints,objective);
 
  


for j = 1:Nr
    
    uq=double(uqtil{j});
    uqref=[uqref,uq];
    
    y=double(ytil{j});
    yref=[yref,y];
    
    
end    

LineVoltageReference = yref+V_tildaProfile;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%% Stage 1 Set up Reference Values for Stage 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uqref=[uqref,uqref(:,1:N)];
yref=[yref,yref(:,1:N)];
xref = zeros(nx,Nr+N+1);
u2ref = zeros(nu2,(Nr+N));
V_tildaProfileExtra = [V_tildaProfile,V_tildaProfile(:,1:N)];

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stage 2 MPC Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R1 = eye(ny);
R2uq = (1/qmax^2)*eye(nuq);
R2u2 = 0*(1/cmax^2)*eye(nu2);
R3 = 100* eye(nx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stage 2 Define Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u2bar = sdpvar(repmat(nu2,1,N),repmat(1,1,N));
uqbar = sdpvar(repmat(nuq,1,N),repmat(1,1,N));
xbar = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
ybar = sdpvar(repmat(ny,1,N),repmat(1,1,N));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stage 2 Setup Optimzation Problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
constraints =[];
objective = 0;
for k = 1:N
    objective = objective + (ybar{k}- yref(:,k))'*R1*(ybar{k}-yref(:,k)) + ...
        (u2bar{k}-u2ref(:,k))'*R2u2*(u2bar{k}-u2ref(:,k)) + ...
        (uqbar{k}-uqref(:,k))'*R2uq*(uqbar{k}-uqref(:,k)) + ...
        (xbar{k}-xref(:,k))'*R3*(xbar{k}-xref(:,k));
    constraints = [constraints, xbar{k+1} == A*xbar{k}+B(numBus+1:end)'.*u2bar{k}];
    constraints = [constraints, ybar{k} == D(:,1:nuq)*uqbar{k}+D(:,(nuq+1):end)*u2bar{k}];
    constraints = [constraints, SOC_des - xbar{k}<= SOC_maxpu];
    constraints = [constraints,  0 <= xbar{k}];
    constraints = [constraints, 0 <= u2bar{k} <= cmax];
    constraints = [constraints, qmin <= uqbar{k} <= qmax];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stage 2 Initialize  Control Desing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% X already taken
xs = [];
U2 = [];
Uq = [];
LineVoltage = [];

xs(:,1) = SOC_initial_conditions;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stage 2 simulate MPC controller
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:Nr
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Stage 2 Solve MPC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    constraintsLoop = constraints;
    for j = 1:N
        constraintsLoop = [constraintsLoop, Vmin*ones(ny,1)-V_tildaProfileExtra(:,k+j-1) <= ybar{j}-Vnom*ones(ny,1) <= Vmax*ones(ny,1)-V_tildaProfileExtra(:,k+j-1)];
    end
    constraintsLoop = [constraintsLoop,xbar{1}==xs(:,k)];
   
    
    diag=optimize(constraintsLoop,objective);
    diagnostics=diag.problem;
    
    % process output
    if diagnostics == 1
        error('The problem is infeasible');
    end
    
    u2=double(u2bar{1});
    U2=[U2,u2];
    
    uq=double(uqbar{1});
    Uq=[Uq,uq];
    
    y=double(ybar{1});
    LineVoltage=[LineVoltage,y+V_tildaProfile(:,k)];
    
    x=double(xbar{2});
    xs=[xs,x];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Full DistFlow Line Voltage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% this section not very important
% RealDistFlowLineVoltage = zeros(ny,Nr); Ynl = zeros(ny,Nr); for i = 1:Nr
%     Qg = Uq(:,i); Pv = U2(:,i); PQv = distFlowSolve(Qg,Pv);
%     RealDistFlowLineVoltage(:,i) = PQv(1:ny,3)/23e3;
% end
 

%% Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
for i = 1:9
    plot(xs(i,:))
    hold on
end
title('xs')
xlabel('timeSteps')
ylabel('SOC_{des} - SOC_{cur}');
legend('EV1','EV2','EV3','EV4','EV5','EV6','EV7','EV8','EV9')
hold off


figure
for i = 1:9
    plot(LineVoltage(i,:));
    hold on
end
title('Line Voltage')
xlabel('timeSteps')
ylabel('Line Voltage (p.u)');
legend('BUS1','BUS2','BUS3','BUS4','BUS5','BUS6','BUS7','BUS8','BUS9')
hold off


% figure
% for i = 1:9
%     plot(RealDistFlowLineVoltage(i,:));
%     hold on
% end
% title('Real DistFlow Line Voltage')
% xlabel('timeSteps')
% ylabel('Nonlinear Model Line Voltage (p.u)');
% hold off

figure
for i = 1:9
    plot(U2(i,:)/1e3)
    hold on
end
title('EV Charging Control')
xlabel('timeSteps')
ylabel('Charging Power in kW');
legend('EV1','EV2','EV3','EV4','EV5','EV6','EV7','EV8','EV9')
hold off

figure
for i = 1:9
    plot(Uq(i,:)/1e3)
    hold on
end
title('Reactive Power Generation Control')
xlabel('timeSteps')
ylabel('Reactive Power Generation in kW');
legend('BUS1','BUS2','BUS3','BUS4','BUS5','BUS6','BUS7','BUS8','BUS9')
hold off

%% Compare with Uncoordinated Charging

% Uncoordinated charging is defined to be charging vehicles at full power
% without any Qg or some constant Qg.
for k = 1:Nr
    y_uncoordinated(:,k) = D(:,(nuq+1):end)*qmax*ones(9,1);
end

LineVoltage_uncoordinated = y_uncoordinated + V_tildaProfile;

figure
for i = 1:9
    plot(LineVoltage_uncoordinated(i,:));
    hold on
end
title('Line Voltage Uncoordinated')
xlabel('timeSteps')
ylabel('Line Voltage (p.u)');
legend('BUS1','BUS2','BUS3','BUS4','BUS5','BUS6','BUS7','BUS8','BUS9')
hold off