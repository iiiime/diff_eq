
function [y,p,i,paras,Tactive,L_seq] = RK_Gradient(paras,y0,h,tend,X0,Y0,g,Edge)

% Simulate for a liear gradient search using Runge-Kutta method (instead of
% Euler)
% Vector y: % Vector y: T1 | T2 | T3 | T4 | T5 | Ap | Yp | Bp | Tau |
% NEXT Tumble? | angles | NEXT Tumble Direction | current x coordinate | current y coordiate | 
% NEXT Mov Direction (0:right;1:down;2:left;3:up) | ChemStablize? |
% BoundaryTumble (override the original tumble?)
% p: probability of receptor activation
% i: iteration step
% paras: parameters used to simulate
% T_active: active form of ligand receptor
% L_seq: concentration of ligand L

% INPUT:
% paras: parameters
% y0: initial condition
% h: step size
% tend: end of the time
% X0: food source-x
% Y0: food source-y
% g: function applied to distance
% Edge: the searching Edge, [xmin xmax ymin ymax]

% When enzyme concentration stabilizes, the bacteria start to move, 
% Kc=steady state Yp
% Stop after the bacteria is with in X0+-3 and Y0+-3

% Add boundary to the searching area


R0 = paras.R0;
% Methylated related constants
Km0 = paras.Km0; % mM
Km1 = paras.Km1; 
Km2 = paras.Km2;
Km3 = paras.Km3;
Km4 = paras.Km4;
Vm0 = paras.Vm0;
Vm1 = paras.Vm1;
Vm2 = paras.Vm2;
Vm3 = paras.Vm3;
Vm4 = paras.Vm4;
% Hill curves
Hm = paras.Hm;
Hc = paras.Hc;
Kc = paras.Kc; % muM
% constant
kR = paras.kR; % /s
kB = paras.kB; %/s
k2B = paras.k2B; % /muM
% kA = 50; % /muM, from appendix
% kY = 100; % /muM, from appendix
kA = paras.kA;
kY = paras.kY;
% Kinetic constants
kZ = paras.kZ; % maybe wrong
KR = paras.KR; % muM
KB = paras.KB; % muM
gammaB = paras.gammaB; % /(muM*s)
gammaY = paras.gammaY;
% total concentration muM
AT = paras.AT;
BT = paras.BT;
% RT = paras.RT;
YT = paras.YT;
ZT = paras.ZT;
% constant cheR
beta = paras.beta;
% tumbling angle distribution
shape = paras.shape;
scale = paras.scale;
% angles = paras.angles; % gamma distribution
% Ligand concentration and activation probability
L = paras.L; % mM



% Moving directions
Dir = {};
Dir{1} = [1;0];
Dir{2} = [0;-1];
Dir{3} = [-1;0];
Dir{4} = [0;1];

% Part of vector y: T1 | T2 | T3 | T4 | T5 | Ap | Yp | Bp | 
f = @(t,y,p) [kB*y(8)*p(2)*y(2)/(KB+p*y(1:5)) - kR*R0*y(1)/(KR+sum(y(1:5)));...
            kR*R0*y(1)/(KR+sum(y(1:5))) + kB*y(8)*p(3)*y(3)/(KB+p*y(1:5))-kR*R0*y(2)/(KR+sum(y(1:5))) - kB*y(8)*p(2)*y(2)/(KB+p*y(1:5));...
            kR*R0*y(2)/(KR+sum(y(1:5))) + kB*y(8)*p(4)*y(4)/(KB+p*y(1:5))-kR*R0*y(3)/(KR+sum(y(1:5))) - kB*y(8)*p(3)*y(3)/(KB+p*y(1:5));...
            kR*R0*y(3)/(KR+sum(y(1:5))) + kB*y(8)*p(5)*y(5)/(KB+p*y(1:5))-kR*R0*y(4)/(KR+sum(y(1:5))) - kB*y(8)*p(4)*y(4)/(KB+p*y(1:5));...
            kR*R0*y(4)/(KR+sum(y(1:5))) - kB*y(8)*p(5)*y(5)/(KB+p*y(1:5));...
            kA*(AT-y(6))*(p*y(1:5))-kY*y(6)*(YT-y(7))-k2B*y(6)*(BT-y(8));...
            kY*y(6)*(YT-y(7)) - kZ*y(7)*ZT - gammaY*y(7);...
            k2B*y(6)*(BT-y(8)) - gammaB*y(8)];

iterations = (tend-h)/h+1;
t = [0:h:(tend-h)];
y = nan(17,size(t,2));
L_seq = nan(1,size(t,2));
Tactive = nan(1,size(t,2));
y(:,1) = y0;
y(13:14,2) = [0;0]; 
StableCheck = 0;

% Vector y: % Vector y: T1 | T2 | T3 | T4 | T5 | Ap | Yp | Bp | Tau |
% NEXT Tumble? | angles | NEXT Tumble Direction | current x coordinate | current y coordiate | 
% NEXT Mov Direction (0:right;1:down;2:left;3:up) | ChemStablize? |
% BoundaryTumble (override the original tumble?)


for i = 1:(iterations-2)
   
    % simulate for enzyme
%     p = [];
%     curr_L = L(abs(y(13,i+1))+1);
  
%     R = 1;
%     D = 1;
%     tau = 2500.0;
%     V = 0;
%     lmd = sqrt(D*tau / (1 + V^2 * tau / (4 * D)));
%     dis = sqrt((y(13,i+1)-x0)^2 + (y(14,i+1)-y0)^2);
%     curr_L = R * exp((y0-y(14,i+1))*V/(2*D))* besseli(0,dis / lmd) / (2 * pi * D);

    dis = sqrt((y(13,i+1)-X0)^2 + (y(14,i+1)-Y0)^2);    
    curr_L = g(dis,y(14,i+1)); % with wind
    % curr_L = g(dis);
    p(1) = Vm0*(1-curr_L^Hm/(curr_L^Hm+Km0^Hm));
    p(2) = Vm1*(1-curr_L^Hm/(curr_L^Hm+Km1^Hm));
    p(3) = Vm2*(1-curr_L^Hm/(curr_L^Hm+Km2^Hm));
    p(4) = Vm3*(1-curr_L^Hm/(curr_L^Hm+Km3^Hm));
    p(5) = Vm4*(1-curr_L^Hm/(curr_L^Hm+Km4^Hm));    
    L_seq(i) = curr_L;
    
    k1 = f(t(i),y(1:8,i),p);
    k2 = f(t(i)+h/2, y(1:8,i)+h/2*k1,p);
    k3 = f(t(i)+h/2, y(1:8,i)+h/2*k2,p);
    k4 = f(t(i)+h, y(1:8,i)+h*k3,p);
    y(1:8,i+1) = y(1:8,i) + h/6*(k1+2*k2+2*k3+k4);   
    Tactive(i+1) = p*y(1:5,i+1);
    if any(y(1:8,i+1)<0) || ...
        y(6,i+1)>AT || y(7,i+1)>YT || y(8,i+1)>BT
        break
    end
    
    % Movement related path    
    y(9,i+1) = y(7,i+1)^paras.Hc / (y(7,i+1)^paras.Hc + paras.Kc^paras.Hc);
    y(10,i+1) = rand(1) < y(9,i+1);
    y(11,i+1) = gamrnd(paras.shape,paras.scale,1,1);
    y(12,i+1) = round(y(11,i+1) / 90);

    if y(10,i+1) == 0
    y(15,i+1) = y(15,i); % remain the origial moving direction
    elseif y(10,i+1) == 1
    y(15,i+1) = mod((y(15,i) + y(12,i+1)),4);
    end
    
    if StableCheck == 0 % later fluctuation will not stp moving
        % after stabilization, then move
        if max(abs(diff(y(7,max(1,i-9):i+1))))>y(7,i)/1000
            y(16,i+1) = 0;
            y(13:14,i+2) = y(13:14,i+1);
            % update Kc
        else
            % the first step Yp is stabilized
            paras.Kc = y(7,i+1);
            y(16,i+1) = 1;
            y(13:14,i+2) = y(13:14,i+1)+Dir{y(15,i+1)+1};
            StableCheck = 1;
        end
    else
        y(16,i+1) = 1;
        NextPos = y(13:14,i+1)+Dir{y(15,i+1)+1};
        for count = 1:100
            if NextPos(1)<Edge(1)||NextPos(1)>Edge(2)||NextPos(2)<Edge(3)||NextPos(2)>Edge(4)
                y(15,i+1) = randi(4)-1;
                y(17,i) = 1;
                NextPos = y(13:14,i+1)+Dir{y(15,i+1)+1};
            else
                y(13:14,i+2) = NextPos;
                y(17,i+1) = 0;
                break
            end
        end
    end
       
    if abs(y(13,i+2)-X0)<2 && abs(y(14,i+1)-Y0)<2
        break
    end

    
end

end



