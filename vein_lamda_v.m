% % Use Euler's method to integrate simple one variable ODE
% % Use program as template for more interesting models

% % Note:  produces slightly different output compared with 
% % example shown in class.  Different values of b, dt, tlast
clear all
close all

% define constant parameter values 
%alpha = 50; % /(muM*h)
%delta_A = 4; % /(muM*h)
%delta_G = 4; % /(muM*h)
%D = 0.4; % /h
%D_E = 0.6; % /h
%G_E = 30; % mM
%beta_H = 50; % muM/h
%gamma_H = 7.5; %/h
%K_H = 7.2; % mM
%n = 7; 
beta_r = 0.14; % /(mM*h)
beta_v = 1.5; % need tuning
gamma_r = 3; % /h
eta = 100; % /mM
K_k = 0.5; % mM
m = 2;
lamda_i = 0;
lamda_p = 2; % /h
theta1 = 0;
theta2 = pi;
omega = 2;
K_kv = 0.15;

% define tunable parameters/initial values
LAMDA_V = [0.2:0.2:0.4];
n = numel(LAMDA_V);

% define time
dt = 0.01 ;  
tlast = 20 ; % h
iterations = round(tlast/dt) ;
time = dt*(1:iterations) ;

% define variables
R_i = zeros(iterations,n) ;
R_p = zeros(iterations,n) ;
R_v = zeros(iterations,n);
RO_i = zeros(iterations,n);
RO_p = zeros(iterations,n);
RO_v = zeros(iterations,n);
MU_i = zeros(iterations,n);
MU_p = zeros(iterations,n);
MU_v = zeros(iterations,n);

% calculate known(prefixed) variables
N_i = 0.25 + 0.125*sin(theta1 + omega*time);
N_p = 1 + 0.5*sin(theta2 + omega*time);
K_i = N_i.^m ./ (K_k.^m + N_i.^m);
K_p = N_p.^m ./ (K_k.^m + N_p.^m);
K_v = N_i.^m ./ (K_kv.^m + N_i.^m);

%define initial values;
R_p0 = [1]; % mM
R_i0 = [0.8] ; % mM 
R_v0 = [0.8];
RO_p0 = [0.2];
RO_i0 = [0.1];
RO_v0 = [0.04];
MU_i0 = [1];
MU_p0 = [1];
MU_v0 = [0];


for ii = 1:n

    % add value to tunable parameter
    lamda_v = LAMDA_V(ii);
    
    % add initial value to the variable vector
    R_i(1,ii) = R_i0(1);
    R_p(1,ii) = R_p0(1);
    R_v(1,ii) = R_v0(1);
    RO_p(1,ii) = RO_p0(1);
    RO_i(1,ii) = RO_i0(1);
    RO_v(1,ii) = RO_v0(1);
    MU_i(1,ii) = MU_i0(1);
    MU_p(1,ii) = MU_p0(1);
    MU_v(1,ii) = MU_v0(1);
    
    for i = 1: (length(time)-1)
        ni = N_i(i);
        np = N_p(i);
        ki = K_i(i);
        kp = K_p(i);
        kv = K_v(i);
        ri = R_i(i,ii);
        rp = R_p(i,ii);
        rv = R_v(i,ii);
        rop = RO_p(i,ii);
        roi = RO_i(i,ii);
        rov = RO_v(i,ii);
        mui = MU_i(i,ii);
        mup = MU_p(i,ii);
        muv = MU_v(i,ii);
        
        dri = beta_r*ni^2 - gamma_r*ri;
        drp = beta_r*np^2 - gamma_r*rp;
        drv = beta_v*ni^2 - gamma_r*rv;
        droi = eta*ri*roi*(1-roi/ki) - lamda_i*roi;
        drop = eta*rp*rop*(1-rop/kp) - lamda_p*rop;
        drov = eta*ri*rov*(1-rov/kv) - lamda_v*rov;
        
        R_i(i+1,ii) = ri + dri*dt;
        R_p(i+1,ii) = rp + drp*dt;
        R_v(i+1,ii) = rv + drv*dt;
        RO_i(i+1,ii) = roi + droi*dt;
        RO_p(i+1,ii) = rop + drop*dt;
        RO_v(i+1,ii) = rov + drov*dt;
        MU_i(i+1,ii) = eta*R_i(i+1,ii)*RO_i(i+1,ii)*(1-RO_i(i+1,ii)/K_i(i+1));
        MU_p(i+1,ii) = eta*R_p(i+1,ii)*RO_p(i+1,ii)*(1-RO_p(i+1,ii)/K_p(i+1));
        MU_v(i+1,ii) = eta*R_v(i+1,ii)*RO_v(i+1,ii)*(1-RO_v(i+1,ii)/K_v(i+1));
        
    end
    ii
end

%% Plotting %%%%%

figure(1)
clf

subplot(2,2,1)
hold on
for i = 1
    plot(time, MU_i(:,i), 'color',[1 0 0,0.25],'linewidth',2)
    plot(time, MU_p(:,i), 'color',[1 0 0],'linewidth',2)
    plot(time, MU_v(:,i),'color',[0 0 0,0.25],'linewidth',2)
end
ylim([-1 5])
ylabel('Growth Rate');
legend({'Interior Growth Rate', 'Peripheral Growth Rate','Vein Cell Growth Rate'})
hold off
set(gca, 'fontsize',15);

subplot(2,2,2)
hold on 
plot(time, K_i, 'color', [0 0 1,0.25], 'linewidth',2)
plot(time,K_p,'color',[0 0 1],'linewidth',2)
plot(time, K_v,'color',[0 0 0,0.25],'linewidth',2)
legend({'Interior', 'Peripheral','Vein Cell'})
ylabel('K(carrying capacity)')
xlabel('Time (h)');
ylim([0 1.5])
hold off
set(gca, 'fontsize',15);

subplot(2,2,3)
hold on
for i = 1
    plot(time, RO_i(:,i), 'color',[1 0 0,0.25],'linewidth',2)
    plot(time, RO_p(:,i), 'color',[1 0 0],'linewidth',2)
    plot(time, RO_v(:,i),'color',[0 0 0,0.25],'linewidth',2)
end
legend({'Interior', 'Peripheral','Vein Cell'})
ylim([0 1.5])
hold off
ylabel('Density')
xlabel('Time(h)')
set(gca, 'fontsize',15);

subplot(2,2,4)
hold on
for i = 1
    plot(time,R_i(:,i),'color',[1 0 0,0.25],'linewidth',2)
    plot(time,R_p(:,i),'color',[1 0 0],'linewidth',2)
    plot(time, R_v(:,i),'color',[0 0 0,0.25],'linewidth',2)
end
legend({'Interior','Peripheral','Vein Cell'})
ylabel('R(ribosomal protein')
xlabel('Time(h)')
ylim([0 0.4])
hold off
set(gca, 'fontsize',15)

%%

figure(2)
hold on
plot(time,RO_p(:,1),'color',[1 0 0,0.25],'linewidth',2)
plot(time,K_p,'color',[0 0 1,0.25],'linewidth',2)
plot(time,RO_v(:,1),'color',[1 0 0],'linewidth',2)
plot(time,K_v,'color',[0 0 1],'linewidth',2)
legend({'Peripheral density','Peripheral K','Vein density','Vein K'})
hold off
