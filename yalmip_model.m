clear
%To use the code, yalmip must be installed:
%In the command window:
%tbxmanager install yalmip
%tbxmanager restorepath
%tbxmanager savepath
yalmip('clear');
% Loading data (loads)
load dati_new.mat

Ng=3; %Number of power plants
Nt=1; %Number of thermic plants
Nl=4; %Number of loads
Nlt=4; %Number of thermic loads
Nls=1; %Number of shiftable loads
% cost function alpha and beta coefficients
betagrid=prezzi_medi(:,2);
alpha1= 0.026; beta1=0.2;      % 
alpha2=0.0446; beta2=2.6;  % 
alpha3=0.018; beta3=0.55;     % 
alpha=diag([alpha1;alpha2;alpha3]);
beta=[beta1;beta2;beta3];

% PV
PV=zeros(24,1);

% Wind
Wind=zeros(24,1);

%Csp_input
Csp_input=Csp_input*1e3+2000; %vedi dati_new.mat

% cost function alpha and beta coefficients
alphaT=0.069;
betaT=3.6;

% L is an array (Nl,24)
L=hotel(:,1)+ospedale(:,1)+ufficio(:,1)+residenziale(:,1);
%L=Potenze_attive_carico.';
% Thermic load
%LT=Potenze_termiche.';
LT=hotel(:,2)+ospedale(:,2)+ufficio(:,2)+residenziale(:,2);
% LS is the array of shiftable loads (Nls,1)
LS = 0;
% Constraints on co-generators power
Gmax=repmat([100;150;80],1,24);
Gmin=repmat([0;-150;0],1,24);
Gridmax=repmat(10000,1,24);
Gridmin=repmat(-10000,1,24);

% Constraints on thermic power to electric grid
THEmax=repmat([5000],1,24);
THEmin=repmat([0],1,24);
% Constraints on thermic power to thermic grid
THTmax=repmat([5000],1,24);
THTmin=repmat([0],1,24);

HH=[0;1;0]; %Capacity storage (0-inactive, 1-active)
HHT=[0]; %Thermic Capacity storage (0-inactive, 1-active)

E0=[0;0;0]; %Initial condition on capacity storage (energy)
Emax=[0;150;0]; %Maximum capacity storage
ETH0=1e8*[1]; %Initial condition on thermic capacity storage (energy)
ETHmax=1e10*[1]; %Maximum thermic capacity storage

%Design variables
act=binvar(Nls,24);
G=sdpvar(Ng,24); %electric power
Grid=sdpvar(1,24); %Grid electric power
E=sdpvar(Ng,24); %capacity storage
THE=sdpvar(1,24); %Csp thermic power
THT=sdpvar(1,24); %Csp thermic power
ETH=sdpvar(Nt,24); %thermic storage
etaT=[1]; %thermic coefficient (thermic -> thermic)
etaE=[0.38]; %thermic coefficient (thermic -> electric)
constraints=[];
%Constraints on shiftable loads:
%act is a matrix of binary design variables {0,1} to indicate the time of
%switching the load
for i=1:Nls
    constraints=[constraints;sum(act(i,:))==1];
end
E(:,1)=E0;
ETH(:,1)=ETH0;
J=0;
prova=sdpvar(1);
for i=1:24
    if i<24
        %Constraint to define the dynamics of energy storage
        constraints=[constraints;E(:,i+1)==E(:,i)-HH.*G(:,i)]; 
        constraints=[constraints;ETH(:,i+1)==ETH(:,i)-HHT.*(THE(:,i)+THT(:,i))]; 
        %Constraints on maximum and minimum capacity storage
        constraints=[constraints;E(:,i+1)>=0;E(:,i+1)<=HH.*Emax];
        constraints=[constraints;HHT.*ETH(:,i+1)>=0;HHT.*ETH(:,i+1)<=HHT.*ETHmax];
    else
        %Constraint to unload the capacity storage at 24th hour
        %constraints=[constraints;E(:,24)==0];
    end
    %Constraint on co-generator power (minimum and maximum)
    constraints=[constraints;G(:,i)<=Gmax(:,i);G(:,i)>=Gmin(:,i)];
    constraints=[constraints;Grid(1,i)<=Gridmax(1,i);Grid(1,i)>=Gridmin(1,i)];
    constraints=[constraints;THT(:,i)<=THTmax(:,i);THT(:,i)>=THTmin(:,i)];
    constraints=[constraints;THE(:,i)<=THEmax(:,i);THE(:,i)>=THEmin(:,i)];
    %Equality constraint
    constraints=[constraints;sum(G(:,i))+Grid(1,i)+sum(etaE.*THE(:,i))+PV(i)+Wind(i)==sum(L(i))+sum(act(:,i).*LS)];
    constraints=[constraints;sum(etaT.*THT(:,i))==LT(i)];%+1e-2;sum(etaT.*THT(:,i))>=LT(i)-1e-2];
    %Csp
    constraints=[constraints;sum(THT(:,i))+sum(THE(:,i))<=Csp_input(i)];
    %Objective function
    J=J+G(:,i).'*alpha*G(:,i)+beta'*G(:,i)+betagrid(i)*Grid(1,i)+...
        (THT(:,i)+THE(:,i)).'*alphaT*(THT(:,i)+THE(:,i))+betaT'*(THT(:,i)+THE(:,i));
end
% options: currently the optimization is based on fmincon, to avoid cplex
options = sdpsettings('solver','fmincon');
% yalmip function to generate an object "optimizer" (prova is only a dummy variable. Do not delete it!).
opt1=optimizer(constraints,J,options,[prova],{G,Grid,E,THT,THE,ETH,J,act});
% Optimizer launch:
% Results are in the variable "res", a cell array
% The following function can be launched by an interpreted matlab function
% in simulink.
[res, errorcode] = opt1(0);
% Gres are the optimized powers
% Eres is the dynamics of energy storage
% Gres are the optimized thermic powers
% Eres is the dynamics of thermic energy storage
% Jres is the objective function
% act_res is the array to indicate the activation of shiftable loads
Gres=res{1};
Gridres=res{2};
Eres=res{3};
THTres=res{4};
THEres=res{5};
ETHres=res{6};
Jres=res{7};
act_res=res{8};
%% for Simulink FromWorkspace blocks
time=3600*[0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 18 19 19 20 20 21 21 22 22 23 23 24];
time1=[1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 18 19 19 20 20 21 21 22 22 23 23 24 24];
Gres_time.time=time;
Gres_time.signals.values=Gres(:,time1).';
Gres_time.signals.dimensions=size(Gres,1);
Gridres_time.time=time;
Gridres_time.signals.values=Gridres(:,time1).';
Gridres_time.signals.dimensions=size(Gridres,1);
Eres_time.time=time;
Eres_time.signals.values=Eres(:,time1).';
Eres_time.signals.dimensions=size(Eres,1);
THTres_time.time=time;
THTres_time.signals.values=THTres(:,time1).';
THTres_time.signals.dimensions=size(THTres,1);
THEres_time.time=time;
THEres_time.signals.values=THEres(:,time1).';
THEres_time.signals.dimensions=size(THEres,1);
ETHres_time.time=time;
ETHres_time.signals.values=ETHres(:,time1).';
ETHres_time.signals.dimensions=size(ETHres,1);
