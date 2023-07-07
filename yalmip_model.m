clear
%To use the code, yalmip must be installed:
%In the command window:
%tbxmanager install yalmip
%tbxmanager restorepath
%tbxmanager savepath
yalmip('clear');
Ng=3; %Number of power plants (co-generator power plant - csp)
Nl=1; %Number of loads
Nls=1; %Number of shiftable loads
% cost function alpha and beta coefficients
alpha1= 0.0021/1000; beta1=34.4/1000;      % 
alpha2=(0.042/1000); beta2= (25.78/1000);  % 
alpha3= 0; beta3=1e-4;     % 
alpha=diag([alpha1;alpha2;alpha3]);
beta=[beta1;beta2;beta3];

% Loading data (loads)
% L is an array (Nl,24)
load dati.mat
L=Potenze_attive_carico.';
% LS is the array of shiftable loads (Nls,1)
LS=200;
% Constraints on co-generators power
Gmax=repmat([10000;100;1000],1,24);
Gmin=repmat([-10000;100;-1000],1,24);

HH=[0;0;1]; %Capacity storage (0-inactive, 1-active)

E0=[0;0;0]; %Initial condition on capacity storage (energy)
Emax=[0;0;10000]; %Maximum capacity storage

%Design variables
act=binvar(Nls,24);
G=sdpvar(Ng,24); %electric power
E=sdpvar(Ng,24); %capacity storage

constraints=[];
%Constraints on shiftable loads:
%act is a matrix of binary design variables {0,1} to indicate the time of
%switching the load
for i=1:Nls
    constraints=[constraints;sum(act(i,:))==1];
end
E(:,1)=E0;
J=0;
prova=sdpvar(1);
for i=1:24
    if i<24
        %Constraint to define the dynamics of energy storage
        constraints=[constraints;E(:,i+1)==E(:,i)-HH.*G(:,i)]; 
        %Constraints on maximum and minimum capacity storage
        constraints=[constraints;HH.*E(:,i+1)>=0;HH.*E(:,i+1)<=Emax];
    else
        %Constraint to unload the capacity storage at 24th hour
        constraints=[constraints;E(:,24)==0];
    end
    %Constraint on co-generator power (minimum and maximum)
    constraints=[constraints;G(:,i)<=Gmax(:,i);G(:,i)>=Gmin(:,i)];
    %Equality constraint
    constraints=[constraints;sum(G(:,i))==sum(L(:,i))+sum(act(:,i).*LS)];
    %Objective function
    J=J+G(:,i).'*alpha*G(:,i)+beta'*G(:,i);
end
% options: currently the optimization is based on fmincon, to avoid cplex
options = sdpsettings('solver','fmincon');
% yalmip function to generate an object "optimizer" (prova is only a dummy variable. Do not delete it!).
opt1=optimizer(constraints,J,options,[prova],{G,E,J,act});
% Optimizer launch:
% Results are in the variable "res", a cell array
% The following function can be launched by an interpreted matlab function
% in simulink.
[res, errorcode] = opt1(0);
% Gres are the optimized powers
% Eres is the dynamics of energy storage
% Jres is the objective function
% act_res is the array to indicate the activation of shiftable loads
Gres=res{1};
Eres=res{2};
Jres=res{3};
act_res=res{4};