clear
%To use the code, yalmip must be installed:
%In the command window:
%tbxmanager install yalmip
%tbxmanager restorepath
%tbxmanager savepath
yalmip('clear');
% Loading data (loads)
load dati.mat

Ng=3; %Number of power plants
Nt=3; %Number of thermic plants
Nl=1; %Number of loads
Nlt=1; %Number of thermic loads
Nls=1; %Number of shiftable loads
% cost function alpha and beta coefficients
betagrid=prezzi_medi;
alpha1= 0.0021/1000; beta1=34.4/1000;      % 
alpha2=(0.042/1000); beta2= (25.78/1000);  % 
alpha3= 0; beta3=1e-4;     % 
alpha=diag([alpha1;alpha2;alpha3]);
beta=[beta1;beta2;beta3];

% cost function alpha and beta coefficients
alphaT=diag([alpha1;alpha2;alpha3]);
betaT=[beta1;beta2;beta3];

% L is an array (Nl,24)
L=Potenze_attive_carico.';
% Thermic load
LT=Potenze_termiche.';
% LS is the array of shiftable loads (Nls,1)
LS =200;
% Constraints on co-generators power
Gmax=repmat([10000;100;1000],1,24);
Gmin=repmat([0;100;-1000],1,24);
Gridmax=repmat(10000,1,24);
Gridmin=repmat(-10000,1,24);

% Constraints on thermic power to electric grid
THEmax=repmat([10000;100;1000],1,24);
THEmin=repmat([0;0;0],1,24);
% Constraints on thermic power to thermic grid
THTmax=repmat([10000;100;1000],1,24);
THTmin=repmat([0;0;0],1,24);

HH=[0;0;1]; %Capacity storage (0-inactive, 1-active)
HHT=[1;1;1]; %Thermic Capacity storage (0-inactive, 1-active)

E0=[0;0;0]; %Initial condition on capacity storage (energy)
Emax=[0;0;10000]; %Maximum capacity storage
ETH0=1e8*[1;1;1]; %Initial condition on thermic capacity storage (energy)
ETHmax=1e8*[1;1;1]; %Maximum thermic capacity storage

%Design variables
act=binvar(Nls,24);
G=sdpvar(Ng,24); %electric power
Grid=sdpvar(1,24); %Grid electric power
E=sdpvar(Ng,24); %capacity storage
THE=sdpvar(Ng,24); %thermic power
THT=sdpvar(Ng,24); %thermic power
ETH=sdpvar(Ng,24); %thermic storage
etaT=[1;0;0]; %thermic coefficient (thermic -> thermic)
etaE=[0;0;0.9]; %thermic coefficient (thermic -> electric)
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
        constraints=[constraints;HH.*E(:,i+1)>=0;HH.*E(:,i+1)<=Emax];
        constraints=[constraints;HHT.*ETH(:,i+1)>=0;HHT.*ETH(:,i+1)<=ETHmax];
    else
        %Constraint to unload the capacity storage at 24th hour
        constraints=[constraints;E(:,24)==0];
    end
    %Constraint on co-generator power (minimum and maximum)
    constraints=[constraints;G(:,i)<=Gmax(:,i);G(:,i)>=Gmin(:,i)];
    constraints=[constraints;Grid(1,i)<=Gridmax(1,i);Grid(1,i)>=Gridmin(1,i)];
    constraints=[constraints;THT(:,i)<=THTmax(:,i);THT(:,i)>=THTmin(:,i)];
    constraints=[constraints;THE(:,i)<=THEmax(:,i);THE(:,i)>=THEmin(:,i)];
    %Equality constraint
    constraints=[constraints;sum(G(:,i))+Grid(1,i)+sum(etaE.*THE(:,i))==sum(L(:,i))+sum(act(:,i).*LS)];
    constraints=[constraints;sum(etaT.*THT(:,i))==sum(LT(:,i))];
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