function dP = toggle_switch2(t,P)

% Parameters
V = 5; % Max translation rate
Kd=10; % Dissociation constant for repressor
n1=4; % Cooperativity of repressor

% Components of ODEs
mRNA1=P(1);
Repressor2=P(2);
mRNA2=P(3);
Repressor1=P(4);
K1=V/(1+(Repressor1/(Kd))^n1); 
K2=V/(1+(Repressor2/(Kd))^n1);

% ODEs
dP(1)=K1-mRNA1; 
dP(2)=V*mRNA1-Repressor2; 
dP(3)=K2-mRNA2; 
dP(4)=V*mRNA2-Repressor1;

dP=dP'; 


end