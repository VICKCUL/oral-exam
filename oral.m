%% HH Model for Oral Exam
clear
close all

T=200;
dt=.001;
time=0:dt:T;

% % Pulsatile input:
t01=50; % Time of pulse in ms
t02 = 100;
t03 = 150;
sigma1=1;% Width of pulse in ms
strength1=6; % Height of pulse
sigma2=1;% Width of pulse in ms
strength2=12; % Height of pulse
sigma3=1;% Width of pulse in ms
strength3=3; % Height of pulse
Iapp1=strength1*exp(-(time-t01).^2/(2*sigma1^2));% gaussian external input
Iapp2=strength2*exp(-(time-t02).^2/(2*sigma2^2));
Iapp3=strength3*exp(-(time-t03).^2/(2*sigma3^2));

Iapp = Iapp1+Iapp2+Iapp3;
% Pulsatile input:
% No current
% Iapp=zeros(size(time));
% t0=150; % Time of pulse in ms
% sigma=1;% Width of pulse in ms
% strength=50; % Height of pulse
% Iapp=Iapp+strength*exp(-(time-t0).^2/(2*sigma^2));

% n variable
alphan=@(V)(.01*(V+55)./(1-exp(-.1*(V+55))));
betan=@(V)(.125*exp(-.0125*(V+65)));
ninfty=@(V)(alphan(V)./(alphan(V)+betan(V)));
taun=@(V)(1./(alphan(V)+betan(V)));

% m variable
alpham=@(V)(.1*(V+40)./(1-exp(-.1*(V+40))));
betam=@(V)(4*exp(-.0556*(V+65)));
minfty=@(V)(alpham(V)./(alpham(V)+betam(V)));
taum=@(V)(1./(alpham(V)+betam(V)));

% h variable
alphah=@(V)(.07*exp(-.05*(V+65)));
betah=@(V)(1./(1+exp(-.1*(V+35))));
hinfty=@(V)(alphah(V)./(alphah(V)+betah(V)));
tauh=@(V)(1./(alphah(V)+betah(V)));

% Parameters
Cm=1;
gL=.3;
EL=-54.387;
gK=36;
EK=-77;
gNa=120;
ENa=50;

% Initial conditions near their fixed points
n0=0.3177;    
m0=0.0530;
h0=0.5960;
V0=-65;

% Currents
IL=@(V)(-gL*(V-EL));
IK=@(n,V)(-gK*n.^4.*(V-EK));
INa=@(m,h,V)(-gNa*m.^3.*h.*(V-ENa));

% Toal ion currents
Iion=@(n,m,h,V)(IL(V)+IK(n,V)+INa(m,h,V));

% Euler solver
V=zeros(size(time));
n=zeros(size(time));
m=zeros(size(time));
h=zeros(size(time));
V(1)=V0;
n(1)=n0;
m(1)=m0;
h(1)=h0;

for i=2:numel(time)
    
    % Update gating variables
    n(i)=n(i-1)+dt*((1-n(i-1)).*alphan(V(i-1))-n(i-1).*betan(V(i-1)));
    m(i)=m(i-1)+dt*((1-m(i-1)).*alpham(V(i-1))-m(i-1).*betam(V(i-1)));
    h(i)=h(i-1)+dt*((1-h(i-1)).*alphah(V(i-1))-h(i-1).*betah(V(i-1)));
    
    % Update membrane potential
    V(i)=V(i-1)+dt*(Iion(n(i-1),m(i-1),h(i-1),V(i-1))+Iapp(i-1))/Cm;

    
end

figure

Vs=-100:.2:0;
subplot(3,1,1)
plot(time,Iapp,'k','linewidth',2)
ylabel('Iext')
set(gca,'LineWidth',1.5)
set(gca,'FontSize',16)
xlabel('time (ms)')
title('inject different Normal strenghs at three locations');

subplot(3,1,2)
plot(time,n,'linewidth',2)
hold on    
plot(time,m,'linewidth',2)
plot(time,h,'linewidth',2)
leg=legend('n','m','h','location','east');
set(leg,'FontSize',12)
ylabel('active probability')
box off
temp=axis;
temp(3)=0;
axis(temp)
set(gca,'LineWidth',1.5)
set(gca,'FontSize',16)
xlabel('time (ms)')
title('Iron Dynamics')

subplot(3,1,3)
plot(time,V,'Color',[.5 .5 .5],'linewidth',2)
set(gca,'LineWidth',1.5)
ylabel('V (mV)')
xlabel('time (ms)')
axis([0 T min(-70,min(V)) max(-57,max(V))])
set(gca,'FontSize',16)
title('Spike Dynamics')

%%





%% LIF model
clear
close all

% Discretized time
T=1000;
dt=.1;
time=0:dt:T;

% Applied current
I=zeros(size(time));
I(time>=50 & time<=200)=10;
I(time>=300 & time<=400)=5;
I(time>=500 & time<=505)=100;
I(time>=600 & time<=800)=-10;
I(time>=900 & time<=905)=100

% Neuron parameters
EL=-72;
taum=15;
Vth = -55;
Vre = -75;


% Euler's method solution
V=zeros(size(time));
rho=zeros(size(time));% keep tracking of the spikes
V(1)=-72;

for i=2:numel(time)
    V(i)=V(i-1)+(1/taum)*(-(V(i-1)-EL)+I(i))*dt;
    if(V(i)>=Vth)
        V(i) = Vre;
        rho(i) = 1/dt;
    end
end

% Plot I and both solutions
figure
subplot(3,1,1)
plot(time,I,'LineWidth',2)
ylabel('I (mV)')
set(gca,'FontSize',18)
title('LIF Model Simulation')

subplot(3,1,2)
plot(time,V,'LineWidth',2)
hold on
ylabel('V (mV)')
set(gca,'FontSize',18)

subplot(3,1,3)
plot(time, rho, 'LineWidth',2)
hold on
ylabel('spike train')
set(gca,'FontSize',18)
xlabel('time (ms)')

%% EIF model
clear
close all

% Discretized time
T=1000;
dt=.1;
time=0:dt:T;

% Applied current
I=13.5+zeros(size(time)); % get 5 spikes
% I=5+zeros(size(time));

% Neuron parameters
EL=-72;
taum=15;
Vth=10;
Vre=-75;
VT=-55;
DT=4;


% Euler's method solution
V=zeros(size(time));
rho=zeros(size(time));
V(1)=-75;
for i=2:numel(time)
    V(i)=V(i-1)+(1/taum)*(-(V(i-1)-EL)+DT*exp((V(i-1)-VT)/DT)+I(i))*dt;
    if(V(i)>=Vth)
        V(i)=Vre;
        V(i-1)=Vth;
        rho(i)=1/dt;
    end
end

% Plot Iapp and both solutions
figure
subplot(3,1,1)
plot(time,I,'LineWidth',2)
ylabel('I (mV)')
set(gca,'FontSize',18)
title('EIF Model Simulation with constant input')

subplot(3,1,2)
plot(time,V,'LineWidth',2)
ylabel('V (mV)')
set(gca,'FontSize',18)

subplot(3,1,3)
plot(time,rho,'LineWidth',2)
ylabel('spike train')
xlabel('time (ms)')
set(gca,'FontSize',18)


%% Stability analysis, strong input derive V to unstable equilibrum, hence spike
clear
close all

% Neuron parameters
EL=-72;
taum=15;
Vth=10;
Vre=-75;
VT=-55;
DT=4;

% Time-constant input strength
% I0=0;
% I0=5;
I0=13.5; % the larger I0, then fixed point disapear, it became saddle point bification.
% RHS of ODE
f=@(V)((-(V-EL)+DT*exp((V-VT)/DT)+I0)/taum);

figure
Vplot=-80:.1:-40;
plot(Vplot,f(Vplot),'LineWidth',2)
hold on
plot(Vplot,0*Vplot,'r--','LineWidth',2)
xlabel('V')
ylabel('f(V)')
set(gca,'FontSize',14)


%% Synaptic EIF model with 1 exci 
% Simulation of an exponential integrate-and-fire (EIF) model neuron
% with current-based snaptic input

clear
close all

% Neuron parameters
EL=-72;
taum=10;
Vth=10;
Vre=-72;
VT=-55;
DT=2;

% Synaptic parameters
taue=4;
taui=6;
Je=80; % here is the differece between conductance based model
Ji=-80;

% Time vector in ms
T=200;
dt=.1;
time=0:dt:T;

% Applied current
Iapp=0+zeros(size(time)); % comparing 

% Excitatory and inhibitory spike trains
te=[20 30 40];
Se=hist(te,time)/dt; % sum of delta functions

% ti=[150];
% ti = [20 40 150 180];
ti = [];
Si=hist(ti,time)/dt;

% Solve ODE with threshold-reset condition
V=zeros(size(time)); % initialize V
V(1)=EL; % Initial condition
Ie=zeros(size(time)); % initialize ge
Ie(1)=0; % Initial condition
Ii=zeros(size(time)); % initialize gi
Ii(1)=0; % Initial condition
N=0;  % Number of spikes
S=zeros(size(time));
for i=2:numel(time)
   
   % Euler step
   V(i)=V(i-1)+dt*((-(V(i-1)-EL)+DT*exp((V(i-1)-VT)/DT)+Ie(i-1)+Ii(i-1)+Iapp(i))/taum);
   Ie(i)=Ie(i-1)+dt*(-Ie(i-1)+(Je/taue)*Se(i-1))/taue;
   Ii(i)=Ii(i-1)+dt*(-Ii(i-1)+(Ji/taui)*Si(i-1))/taui;
   
   if(V(i)>=Vth)   % If membrane potential is above threshold
       N=N+1;
       S(i)=1/dt;
       V(i)=Vre;   % Reset membrane potential
       V(i-1)=Vth; % This makes the plots look nicer
   end
    
    
end

% 3 graphs, since only need to plot current, not the conductance.

% Plot results    
figure
subplot(3,1,1)
plot(time,Ie,'r','LineWidth',2)
ylabel('I_e')
%title("one excitatory neuron input ")
xlabel('times(ms)')
axis tight 
box off
set(gca,'FontSize',15)

figure 
subplot(3,1,2)
plot(time,Ii,'LineWidth',2)
ylabel('I_i')
%title("one inhibitory neuron input ")
xlabel('times(ms)')
axis tight 
box off
set(gca,'FontSize',15)

figure
subplot(3,1,3)
plot(time,V,'k-','LineWidth',2)
ylabel('V (mV)')
xlabel('time (ms)')
%title("postsynaptic neruon current")
axis tight 
box off
set(gca,'FontSize',15)

%% Synaptic EIF model with 1 exci strong input
% Simulation of an exponential integrate-and-fire (EIF) model neuron
% with current-based snaptic input

clear
close all

% Neuron parameters
EL=-72;
taum=10;
Vth=10;
Vre=-72;
VT=-55;
DT=2;

% Synaptic parameters
taue=4;
taui=6;
Je=80; % here is the differece between conductance based model
Ji=-80;

% Time vector in ms
T=200;
dt=.1;
time=0:dt:T;

% Applied current
Iapp=0+zeros(size(time)); % comparing 

% Excitatory and inhibitory spike trains
te=[20 30 40 50 60 70];
Se=hist(te,time)/dt; % sum of delta functions

% ti=[150];
% ti = [20 40 150 180];
ti = [];
Si=hist(ti,time)/dt;

% Solve ODE with threshold-reset condition
V=zeros(size(time)); % initialize V
V(1)=EL; % Initial condition
Ie=zeros(size(time)); % initialize ge
Ie(1)=0; % Initial condition
Ii=zeros(size(time)); % initialize gi
Ii(1)=0; % Initial condition
N=0;  % Number of spikes
S=zeros(size(time));
for i=2:numel(time)
   
   % Euler step
   V(i)=V(i-1)+dt*((-(V(i-1)-EL)+DT*exp((V(i-1)-VT)/DT)+Ie(i-1)+Ii(i-1)+Iapp(i))/taum);
   Ie(i)=Ie(i-1)+dt*(-Ie(i-1)+(Je/taue)*Se(i-1))/taue;
   Ii(i)=Ii(i-1)+dt*(-Ii(i-1)+(Ji/taui)*Si(i-1))/taui;
   
   if(V(i)>=Vth)   % If membrane potential is above threshold
       N=N+1;
       S(i)=1/dt;
       V(i)=Vre;   % Reset membrane potential
       V(i-1)=Vth; % This makes the plots look nicer
   end
    
    
end

% 3 graphs, since only need to plot current, not the conductance.

% Plot results    
figure
subplot(3,1,3)
plot(time,Ie,'r','LineWidth',2)
ylabel(['$I_{syn}(t)$'], 'Interpreter','latex')
%title(['$' latex(f) '$ for $x$ and $y$ in $[-2\pi,2\pi]$'],'Interpreter','latex')
title("neuron k with strong input ")
xlabel('times(ms)')
axis tight 
box off
set(gca,'FontSize',15)

%% Synaptic EIF model with 1 inh 
% Simulation of an exponential integrate-and-fire (EIF) model neuron
% with current-based snaptic input
clear
close all

% Neuron parameters
EL=-72;
taum=10;
Vth=10;
Vre=-72;
VT=-55;
DT=2;

% Synaptic parameters
taue=4;
taui=6;
Je=80; % here is the differece between conductance based model
Ji=-80;

% Time vector in ms
T=200;
dt=.1;
time=0:dt:T;

% Applied current
Iapp=0+zeros(size(time)); % comparing 

% Excitatory and inhibitory spike trains
te=[]; % adding spike times
Se=hist(te,time)/dt; % sum of delta functions

% ti=[150];
ti = [20 40 60];
Si=hist(ti,time)/dt;

% Solve ODE with threshold-reset condition
V=zeros(size(time)); % initialize V
V(1)=EL; % Initial condition
Ie=zeros(size(time)); % initialize ge
Ie(1)=0; % Initial condition
Ii=zeros(size(time)); % initialize gi
Ii(1)=0; % Initial condition
N=0;  % Number of spikes
S=zeros(size(time));
for i=2:numel(time)
   
   % Euler step
   V(i)=V(i-1)+dt*((-(V(i-1)-EL)+DT*exp((V(i-1)-VT)/DT)+Ie(i-1)+Ii(i-1)+Iapp(i))/taum);
   Ie(i)=Ie(i-1)+dt*(-Ie(i-1)+(Je/taue)*Se(i-1))/taue;
   Ii(i)=Ii(i-1)+dt*(-Ii(i-1)+(Ji/taui)*Si(i-1))/taui;
   
   if(V(i)>=Vth)   % If membrane potential is above threshold
       N=N+1;
       S(i)=1/dt;
       V(i)=Vre;   % Reset membrane potential
       V(i-1)=Vth; % This makes the plots look nicer
   end
    
    
end

% 3 graphs, since only need to plot current, not the conductance.

% Plot results    
figure
subplot(3,1,1)
plot(time,Ie,'r','LineWidth',2)
ylabel('I_e')
%title("one excitatory neuron input ")
xlabel('times(ms)')
axis tight 
box off
set(gca,'FontSize',15)

figure 
subplot(3,1,2)
plot(time,Ii,'LineWidth',2)
ylabel('I_i')
%title("one inhibitory neuron input ")
xlabel('times(ms)')
axis tight 
box off
set(gca,'FontSize',15)

figure
subplot(3,1,3)
plot(time,V,'k-','LineWidth',2)
ylabel('V (mV)')
xlabel('time (ms)')
%title("postsynaptic neruon current")
axis tight 
box off
set(gca,'FontSize',15)
%% 1 exc and 1 inh input to 1 post-syn neuron

% Simulation of an exponential integrate-and-fire (EIF) model neuron
% with current-based snaptic input

clear
close all

% Neuron parameters
EL=-72;
taum=10;
Vth=10;
Vre=-72;
VT=-55;
DT=2;

% Synaptic parameters
taue=4;
taui=6;
Je=80; % here is the differece between conductance based model
Ji=-80;

% Time vector in ms
T=200;
dt=.1;
time=0:dt:T;

% Applied current
Iapp=0+zeros(size(time)); % comparing 
% Iapp=10+zeros(size(time)); 
% Iapp=-8+zeros(size(time));

% Excitatory and inhibitory spike trains
% te=[50]; % spike times
te=[50 60 65 100 175]; % adding spike times
Se=hist(te,time)/dt; % sum of delta functions

% ti=[150];
ti = [20 40 150 180];
Si=hist(ti,time)/dt;



% Solve ODE with threshold-reset condition
V=zeros(size(time)); % initialize V
V(1)=EL; % Initial condition
Ie=zeros(size(time)); % initialize ge
Ie(1)=0; % Initial condition
Ii=zeros(size(time)); % initialize gi
Ii(1)=0; % Initial condition
N=0;  % Number of spikes
S=zeros(size(time));
for i=2:numel(time)
   
   % Euler step
   V(i)=V(i-1)+dt*((-(V(i-1)-EL)+DT*exp((V(i-1)-VT)/DT)+Ie(i-1)+Ii(i-1)+Iapp(i))/taum);
   Ie(i)=Ie(i-1)+dt*(-Ie(i-1)+(Je/taue)*Se(i-1))/taue;
   Ii(i)=Ii(i-1)+dt*(-Ii(i-1)+(Ji/taui)*Si(i-1))/taui;
   
   if(V(i)>=Vth)   % If membrane potential is above threshold
       N=N+1;
       S(i)=1/dt;
       V(i)=Vre;   % Reset membrane potential
       V(i-1)=Vth; % This makes the plots look nicer
   end
    
    
end

% 3 graphs, since only need to plot current, not the conductance.

% Plot results    
figure
subplot(3,1,1)
plot(time,Ie,'r','LineWidth',2)
ylabel('I_e')
%xlabel('time (ms)')
title('1 exc. and 1 inh. input to 1 post-syn. neuron')
axis tight 
box off
set(gca,'FontSize',15)

subplot(3,1,2)
plot(time,Ii,'b','LineWidth',2)
ylabel('I_i')
%xlabel('time (ms)')
axis tight 
box off
set(gca,'FontSize',15)

subplot(3,1,3)
plot(time,V,'k','LineWidth',2)
ylabel('V (mV)')
xlabel('time (ms)')
axis tight 
box off
set(gca,'FontSize',15)

%% EIF with Homogenous Poisson Input
% Simulation of an exponential integrate-and-fire (EIF) model neuron
% with current-based snaptic input

clear
close all

% Neuron parameters
EL=-72;
taum=10;
Vth=10;
Vre=-72;
VT=-55;
DT=2;

% Synaptic parameters
taue=4;
taui=6;
Ne=80;
Ni=20;
re=10/1000+zeros(Ne,1);
ri=15/1000+zeros(Ni,1);
Je=25+zeros(1,Ne);
Ji=-25+zeros(1,Ni);



% Time vector in ms
T=1000;
dt=.1;
time=0:dt:T;

% Applied current
Iapp=0+zeros(size(time));

% Excitatory and inhibitory spike trains
Se=zeros(Ne,numel(time));
Si=zeros(Ni,numel(time));
for k=1:Ne
    Se(k,:)=HomogPoisson(re(k),T,dt);
end
for k=1:Ni
    Si(k,:)=HomogPoisson(ri(k),T,dt);
end


% Solve ODE with threshold-reset condition
V=zeros(size(time)); % initialize V
V(1)=EL; % Initial condition
Ie=zeros(size(time)); % initialize ge
Ie(1)=0; % Initial condition
Ii=zeros(size(time)); % initialize gi
Ii(1)=0; % Initial condition
S=zeros(size(time));
for i=2:numel(time)
   
   % Euler step
   V(i)=V(i-1)+dt*((-(V(i-1)-EL)+DT*exp((V(i-1)-VT)/DT)+Ie(i-1)+Ii(i-1)+Iapp(i))/taum);
   Ie(i)=Ie(i-1)+dt*(-Ie(i-1)+Je*Se(:,i-1))/taue;
   Ii(i)=Ii(i-1)+dt*(-Ii(i-1)+Ji*Si(:,i-1))/taui;
   
   if(V(i)>=Vth)   % If membrane potential is above threshold
       S(i)=1/dt;
       V(i)=Vre;   % Reset membrane potential
       V(i-1)=Vth; % This makes the plots look nicer
   end
    
    
end

Iebar=Je*re;
Iibar=Ji*ri;
Ibar=Iebar+Iibar;

mean(S)

% Plot results    
figure
subplot(2,1,1)
plot(time,Ie,'r','LineWidth',1)
hold on
plot(time,Iebar+0*time,'r--','LineWidth',2)
plot(time,Ii,'b','LineWidth',1)
plot(time,Iibar+0*time,'b--','LineWidth',2)
plot(time,Ie+Ii,'k','LineWidth',1)
plot(time,Ibar+0*time,'k--','LineWidth',2)
title('Synaptic Current Input')
ylabel('I')
%leg=legend('Ii','Ie','Ii','location','east');
%set(leg,'FontSize',12)
%xlabel('time (ms)')
axis tight 
box off
set(gca,'FontSize',15)

figure
subplot(2,1,2)
plot(time,V,'LineWidth',2)
%ylabel('V (mV)')
title('Membrane Potential dynamics for Post-syn Neuron j')
ylabel(['$V_j(mV)$'], 'Interpreter','latex')
xlabel('time (ms)')
axis tight 
box off
set(gca,'FontSize',15)

%% Synaptic EIF model with many neurons input
% Simulation of an exponential integrate-and-fire (EIF) model neuron
% with current-based snaptic input

clear
close all

% Neuron parameters
EL=-72;
taum=10;
Vth=10;
Vre=-72;
VT=-55;
DT=2;

% Synaptic parameters
taue=4;
taui=6;
Ne=80;
Ni=20;
re=10/1000+zeros(Ne,1);
ri=15/1000+zeros(Ni,1);
Je=25+zeros(1,Ne);
Ji=-25+zeros(1,Ni);



% Time vector in ms
T=200;
dt=.1;
time=0:dt:T;

% Applied current
Iapp=0+zeros(size(time));

% Excitatory and inhibitory spike trains
Se=zeros(Ne,numel(time));
Si=zeros(Ni,numel(time));
for k=1:Ne
    Se(k,:)=HomogPoisson(re(k),T,dt);
end
for k=1:Ni
    Si(k,:)=HomogPoisson(ri(k),T,dt);
end


% Solve ODE with threshold-reset condition
V=zeros(size(time)); % initialize V
V(1)=EL; % Initial condition
Ie=zeros(size(time)); % initialize ge
Ie(1)=0; % Initial condition
Ii=zeros(size(time)); % initialize gi
Ii(1)=0; % Initial condition
S=zeros(size(time));
for i=2:numel(time)
   
   % Euler step
   V(i)=V(i-1)+dt*((-(V(i-1)-EL)+DT*exp((V(i-1)-VT)/DT)+Ie(i-1)+Ii(i-1)+Iapp(i))/taum);
   Ie(i)=Ie(i-1)+dt*(-Ie(i-1)+Je*Se(:,i-1))/taue;
   Ii(i)=Ii(i-1)+dt*(-Ii(i-1)+Ji*Si(:,i-1))/taui;
   
   if(V(i)>=Vth)   % If membrane potential is above threshold
       S(i)=1/dt;
       V(i)=Vre;   % Reset membrane potential
       V(i-1)=Vth; % This makes the plots look nicer
   end
    
    
end

Iebar=Je*re;
Iibar=Ji*ri;
Ibar=Iebar+Iibar;

mean(S)

% Plot results    
figure
subplot(2,1,1)
plot(time,Ie,'r','LineWidth',1)
hold on
%plot(time,Iebar+0*time,'r--','LineWidth',2)
plot(time,Ii,'b','LineWidth',1)
%plot(time,Iibar+0*time,'b--','LineWidth',2)
%plot(time,Ie+Ii,'k','LineWidth',1)
%plot(time,Ibar+0*time,'k--','LineWidth',2)
ylabel('Current (mV/ms')
title("Many neurons' input current")
axis tight 
box off
set(gca,'FontSize',15)

figure
subplot(2,1,2)
plot(time,V,'LineWidth',2)
hold on
plot(time,VT+0*time,'r--','LineWidth',2)
%ylabel('V (mV)')
ylabel(['$V_j(mV)$'], 'Interpreter','latex')
xlabel('time (ms)')
axis tight 
box off
set(gca,'FontSize',15)
leg=legend('V','VT','location','east');
set(leg,'FontSize',12)

%% Feedforward Network
% Simulation of an exponential integrate-and-fire (EIF) model neuron
% with current-based snaptic input

clear
close all

% Neuron parameters
EL=-72;
taum=10;
Vth=10;
Vre=-72;
VT=-55;
DT=2;

% Synaptic parameters
taue=4;
taui=6;
Ne=1600;
Ni=400;
N=100;
re=5/1000%10/1000;
ri=15/1000;
pe=.1;
pi=.1;
je=40%20;
ji=-40;

% Build connection matrices
Je=je*binornd(1,pe,N,Ne);
Ji=ji*binornd(1,pi,N,Ni);

% Time vector in ms
T=500;
dt=.1;
time=0:dt:T;

% Applied current
Iapp=0+zeros(N,numel(time));

% Excitatory and inhibitory spike trains
Se=zeros(Ne,numel(time));
Si=zeros(Ni,numel(time));
for k=1:Ne
    Se(k,:)=HomogPoisson(re,T,dt);
end
for k=1:Ni
    Si(k,:)=HomogPoisson(ri,T,dt);
end


% Solve ODE with threshold-reset condition
V=zeros(N,numel(time)); % initialize V
V(:,1)=EL; % Initial condition
Ie=zeros(N,numel(time)); % initialize ge
Ie(:,1)=0; % Initial condition
Ii=zeros(N,numel(time)); % initialize gi
Ii(:,1)=0; % Initial condition
S=zeros(N,numel(time));
for i=2:numel(time)
   
   % Euler step
   V(:,i)=V(:,i-1)+dt*((-(V(:,i-1)-EL)+DT*exp((V(:,i-1)-VT)/DT)+Ie(:,i-1)+Ii(:,i-1)+Iapp(:,i))/taum);
   Ie(:,i)=Ie(:,i-1)+dt*(-Ie(:,i-1)+Je*Se(:,i-1))/taue;
   Ii(:,i)=Ii(:,i-1)+dt*(-Ii(:,i-1)+Ji*Si(:,i-1))/taui;
   
   Ispike=find(V(:,i)>=Vth);
   S(Ispike,i)=1/dt;
   V(Ispike,i)=Vre;   % Reset membrane potential
   V(Ispike,i-1)=Vth; % This makes the plots look nicer    
    
end

MeanRate=1000*mean(S(:))

% Mean-field inputs
Iebar=Ne*pe*je*re;
Iibar=Ni*pi*ji*ri;
Ibar=Iebar+Iibar;

% Plot results    
figure

subplot(3,1,1)
plot(time,mean(Ie),'r','LineWidth',1)
hold on
plot(time,Iebar+0*time,'r--','LineWidth',1)
plot(time,mean(Ii),'b','LineWidth',1)
plot(time,Iibar+0*time,'b--','LineWidth',1)
plot(time,mean(Ie+Ii),'k','LineWidth',1)
plot(time,Ibar+0*time,'k--','LineWidth',1)
ylabel('I')
title("Mean Synaptic Current")
axis tight 
box off
set(gca,'FontSize',15)

subplot(3,1,2)
[NeuronIndices,TimeIndices]=find(S);
SpikeTimes=TimeIndices*dt;
plot(SpikeTimes,NeuronIndices,'k.','MarkerSize',5)
hold on
NeuronOneTime = SpikeTimes(NeuronIndices ==1)
NeuronOneSpike = NeuronIndices(NeuronIndices ==1)
plot(NeuronOneTime,NeuronOneSpike,'b.','MarkerSize',15)
NeuronTwoTime = SpikeTimes(NeuronIndices ==50)
NeuronTwoSpike = NeuronIndices(NeuronIndices ==50)
plot(NeuronTwoTime,NeuronTwoSpike,'r.','MarkerSize',15)
xlabel('time (ms)')
ylabel('neuron index')
title("Post-syn neruons spikes in the network")
axis tight 
set(gca,'FontSize',15)

subplot(3,1,3)
plot(time,V(1,:),"b",'LineWidth',2)
hold on
plot(time,V(50,:),"r",'LineWidth',2)
ylabel('V (mV)')
title("Membrane Potential for Post-syn neuron 1 & neuron 50")
axis tight 
box off
set(gca,'FontSize',15)

%% F-I curve
% Simulation of an exponential integrate-and-fire (EIF) model neuron
% with current-based snaptic input

clear
close all


reVals=(3:1:15)/1000;
numVals=numel(reVals);
Ibar=zeros(numVals,1);
r=zeros(numVals,1);
for j=1:numVals
    [j numVals]
    
    % Neuron parameters
    EL=-72;
    taum=10;
    Vth=10;
    Vre=-72;
    VT=-55;
    DT=2;

    % Synaptic parameters
    taue=4;
    taui=6;
    Ne=1600;
    Ni=400;
    N=100;
    re=reVals(j);%10/1000;
    ri=15/1000;
    pe=.1;
    pi=.1;
    je=20;
    ji=-40;

    % Build connection matrices
    Je=je*binornd(1,pe,N,Ne);
    Ji=ji*binornd(1,pi,N,Ni);

    % Time vector in ms
    T=1000;
    dt=.1;
    time=0:dt:T;

    % Applied current
    Iapp=0+zeros(N,numel(time));

    % Excitatory and inhibitory spike trains
    Se=zeros(Ne,numel(time));
    Si=zeros(Ni,numel(time));
    for k=1:Ne
        Se(k,:)=HomogPoisson(re,T,dt);
    end
    for k=1:Ni
        Si(k,:)=HomogPoisson(ri,T,dt);
    end


    % Solve ODE with threshold-reset condition
    V=zeros(N,numel(time)); % initialize V
    V(:,1)=EL; % Initial condition
    Ie=zeros(N,numel(time)); % initialize ge
    Ie(:,1)=0; % Initial condition
    Ii=zeros(N,numel(time)); % initialize gi
    Ii(:,1)=0; % Initial condition
    S=zeros(N,numel(time));
    for i=2:numel(time)

       % Euler step
       V(:,i)=V(:,i-1)+dt*((-(V(:,i-1)-EL)+DT*exp((V(:,i-1)-VT)/DT)+Ie(:,i-1)+Ii(:,i-1)+Iapp(:,i))/taum);
       Ie(:,i)=Ie(:,i-1)+dt*(-Ie(:,i-1)+Je*Se(:,i-1))/taue;
       Ii(:,i)=Ii(:,i-1)+dt*(-Ii(:,i-1)+Ji*Si(:,i-1))/taui;

       Ispike=find(V(:,i)>=Vth);
       S(Ispike,i)=1/dt;
       V(Ispike,i)=Vre;   % Reset membrane potential
       V(Ispike,i-1)=Vth; % This makes the plots look nicer    

    end

    r(j)=1000*mean(S(:));

    % Mean-field inputs
    Iebar=Ne*pe*je*re;
    Iibar=Ni*pi*ji*ri;
    
    Ibar(j)=Iebar+Iibar;


end





plot(Ibar,r,'LineWidth',2)
hold on
xlabel('Ibar')
ylabel('r (Hz)')
set(gca,'FontSize',16)

%% Recurrent Network

clear 
close all

% Time to run simulation and time bin size
% All time should be interepreted as msec
T=250;
dt=.1;

% Number of exc and inh neurons
Ne=4000;
Ni=1000;

% Connection strengths
jee=10;
jei=-40;
jie=30;
jii=-60;

% Connection probabilities
pee=.05;
pei=.05;
pie=.05;
pii=.05;

% Mean-field connectivity
wee=jee*pee*Ne;
wei=jei*pei*Ni;
wie=jie*pie*Ne;
wii=jii*pii*Ni;
W=[wee wei; wie wii];
lambda=eig(W);


% EIF neuron parameters
taum=10;
EL=-72;
Vreset=-72;
VT=-55;
DeltaT=2;
Vth=0;


% Synaptic time constants
taue=8;
taui=4;

% EIF nonlinearity. Defining this function
% makes the Euler code below more concise
psi=@(V)(DeltaT*exp((V-VT)/DeltaT));

% Time vector and number of time bins
time=0:dt:T;
Nt=numel(time);

% Injected currents, defined as functions of time
Iappe=30+0*normcdf(time,500,10);
Iappi=15+0*normcdf(time,500,10);

% Allocate
Ve=zeros(Ne,Nt);
Vi=zeros(Ni,Nt);
Iee=zeros(Ne,Nt);
Iei=zeros(Ne,Nt);
Iie=zeros(Ni,Nt);
Iii=zeros(Ni,Nt);
Se=zeros(Ne,Nt);
Si=zeros(Ni,Nt);

% Random initial condition
% between EL and VT
Ve(:,1)=randn(Ne,1)*5-60;%(VT-EL)/3+EL;%rand(Ne,1)*(VT-EL)+EL;
Vi(:,1)=randn(Ni,1)*5-60;%(VT-EL)/3+EL;%rand(Ni,1)*(VT-EL)+EL;


% Connection matrices
Jee=jee*sparse(binornd(1,pee,Ne,Ne));
Jei=jei*sparse(binornd(1,pei,Ne,Ni));
Jie=jie*sparse(binornd(1,pie,Ni,Ne));
Jii=jii*sparse(binornd(1,pii,Ni,Ni));

tic
fprintf('\nPercent sim complete:\n')
for ii=2:Nt
    % Euler steps
    Ve(:,ii)=Ve(:,ii-1)+dt*(-(Ve(:,ii-1)-EL)+Iee(:,ii-1)+Iei(:,ii-1)+psi(Ve(:,ii-1))+Iappe(ii))/taum;
    Vi(:,ii)=Vi(:,ii-1)+dt*(-(Vi(:,ii-1)-EL)+Iie(:,ii-1)+Iii(:,ii-1)+psi(Vi(:,ii-1))+Iappi(ii))/taum;
    Iee(:,ii)=Iee(:,ii-1)+dt*(-Iee(:,ii-1)+Jee*Se(:,ii-1))/taue;
    Iie(:,ii)=Iie(:,ii-1)+dt*(-Iie(:,ii-1)+Jie*Se(:,ii-1))/taue;
    Iei(:,ii)=Iei(:,ii-1)+dt*(-Iei(:,ii-1)+Jei*Si(:,ii-1))/taui;
    Iii(:,ii)=Iii(:,ii-1)+dt*(-Iii(:,ii-1)+Jii*Si(:,ii-1))/taui;    
    
    % Find where spikes occur, store them, reset potentials
    temp=find(Ve(:,ii)>=Vth); % Find e spikes    
    Se(temp,ii)=1/dt; % Store spikes
    Ve(temp,ii)=Vreset; % Reset Ve    
    Ve(temp,ii-1)=Vth; % This makes the traces look better
    
    temp=find(Vi(:,ii)>=Vth); % Find i spikes    
    Si(temp,ii)=1/dt; % Store spikes
    Vi(temp,ii)=Vreset; % Reset Vi    
    Vi(temp,ii-1)=Vth; % This makes the traces look better

    % Display percent finished every 10%
    if(mod(100*ii/(Nt-1),10)==0)
        fprintf('%d\n',round(100*ii/Nt))
        drawnow;
    end
end
toc
fprintf('\n')

Erate=1000*mean(mean(Se(:,:)))
Irate=1000*mean(mean(Si(:,:)))

% Find exc spike times
[I,J]=find(Se);
eSpikeTimes=J*dt;
eSpikeIndices=I;

% Find inh spike times
[I,J]=find(Si);
iSpikeTimes=J*dt;
iSpikeIndices=I;

%Make raster plots
figure
%subplot(2,1,1)
plot(eSpikeTimes,eSpikeIndices,'r.')
hold on
plot(iSpikeTimes,iSpikeIndices+4000,'b.')
title('Recurrent Network Spiking Activities')
xlabel('time (ms)')
ylabel('neuron index')
%axis([0 T 1 N])
set(gca,'XLim',[50 T])
set(gca,'FontSize',16)


% Compute smoothed rates of exc and inh populations

% Width of kernel
sigma=2;

% Discretization of kernel time
tau=-3*sigma:dt:3*sigma;

% Gaussian kernel
k=normpdf(tau,0,sigma);
k=k/(sum(k*dt));

% Smoothed firing rate
reSmooth=1000*conv(mean(Se),k,'same')*dt;
riSmooth=1000*conv(mean(Si),k,'same')*dt;

% Ignore first 50 ms
TBurn=50;

figure
subplot(2,1,1)
plot(time(time>=TBurn),mean(Iee(:,time>=TBurn))+Iappe(time>=TBurn),'r','LineWidth',2)
hold on
plot(time(time>=TBurn),mean(Iei(:,time>=TBurn)),'b','LineWidth',2)
%plot(time,Iappe(time))
plot(time(time>=TBurn),mean(Iee(:,time>=TBurn))+mean(Iei(:,time>=TBurn))+Iappe(time>=TBurn),'k','LineWidth',2)
legend('Iee+Iappe','Iei','Itotal')
title('Recurrent Network Simulations')
xlabel('time (ms)')
ylabel('currents')
set(gca,'FontSize',16)

% subplot(3,1,2)
% plot(time(time>=TBurn),mean(Iii(:,time>=TBurn))+Iappi(time>=TBurn),'LineWidth',2)
% hold on
% plot(time(time>=TBurn),mean(Iie(:,time>=TBurn)),'LineWidth',2)
% %plot(time,Iappe(time))
% plot(time(time>=TBurn),mean(Iee(:,time>=TBurn))+mean(Iie(:,time>=TBurn))+Iappi(time>=TBurn),'k','LineWidth',2)
% legend('Iii+Iappi','Iie','Itotal')
% title('Recurrent Network Simulations')
% xlabel('time (ms)')
% ylabel('currents')
% set(gca,'FontSize',16)

subplot(2,1,2)
plot(time(time>=TBurn),reSmooth(time>=TBurn),'r','LineWidth',2)
hold on
plot(time(time>=TBurn),riSmooth(time>=TBurn),'b','LineWidth',2)
plot(time(time>=TBurn),Erate + 0*time(time>=TBurn),'r--', 'linewidth',2)
plot(time(time>=TBurn),Irate + 0*time(time>=TBurn),'b--', 'linewidth',2)
xlabel('time (ms)')
ylabel('rate (Hz)')
%set(gca,'YLim',[0 Inf])
legend('e','i')
set(gca,'FontSize',16)


%% f-I curve
clear
close all
 
% Discretized time
T=10000; %10s
dt=.1;
time=0:dt:T;
 
% Neuron parameters
EL=-72;
taum=15;
Vth=10;
Vre=-75;
VT=-55;
DT=4;

I0i = 10:0.1:20;
rate=zeros(size(I0i));

for j=1:numel(I0i)
    I=I0i(j)+zeros(size(time));
    
    V=zeros(size(time));
    rho=zeros(size(time));
    V(1)=-75;
    for i= 2:numel(time)
        V(i)=V(i-1)+(1/taum)*(-(V(i-1)-EL)+DT*exp((V(i-1)-VT)/DT)+I(i))*dt;
        if(V(i)>=Vth)
            V(i)=Vre;
            V(i-1)=Vth;
            rho(i)=1/dt;
        end
    end
    rate(j) =1000 * mean(rho); % in Hz
end

% non-linear fit
Ith = 12.99;
c = 7.759;
f = @(x)(c.*(x-Ith).^(0.6).*(x>Ith)); 

tau = .1;
ratemodel=zeros(size(I0i,2), numel(time));
for j=1:numel(I0i)
    [j numel(I0i)]
    for i = 1:numel(time)
        ratemodel(j,i) = ratemodel(j,i)+dt*(1/tau)*f(I0i(j));
    end
end


x = I0i;
y = rate;
rateplot = mean(ratemodel,2);

figure
plot(x,y,'k','linewidth',2)
hold on
plot(x,f(x),'b','linewidth',2)
plot(I0i, rateplot,'r--','linewidth',1)
title('f-I curve fitted from Spiking & Rate Models')
xlabel('I (mV)')
ylabel('r (Hz)')
legend('Spiking model', 'fit','Rate model')
set(gca,'FontSize',16)
%%
clear
close all

% Number of neurons in each population
Ne=4000;
Ni=1000;
N=Ne+Ni;

% Number of neurons in ffwd layer
%%%Nx=1000;

% Recurrent net connection probabilities
P=[0.05 0.05; 0.05 0.05];

% Ffwd connection probs
%%%Px=[.05; .05];

% Mean connection strengths between each cell type pair
Jm=15*[50 -300; 225 -500]/sqrt(N);
%%%Jxm=15*[180; 135]/sqrt(N);

% Time (in ms) for sim
T=20000;

% Time discretization
dt=.1;

% FFwd spike train rate (in kHz)
%%%rx=15/1000;

% Number of time bins
Nt=round(T/dt);
time=dt:dt:T;

% Synaptic timescales
taux=8;
taue=6;
taui=4;

%%% added the +.4*sqrt(N) and +.3*sqrt(N)
%%% and changed "Iapp" to "X"
% Time constant ffwd input
Xe=zeros(Ne,1)+.4*sqrt(N);
Xi=zeros(Ni,1)+.3*sqrt(N);
%%%

% Generate connectivity matrices
tic
Jee=Jm(1,1)*binornd(1,P(1,1),Ne,Ne);
Jei=Jm(1,2)*binornd(1,P(1,2),Ne,Ni);
Jie=Jm(2,1)*binornd(1,P(2,1),Ni,Ne);
Jii=Jm(2,2)*binornd(1,P(2,2),Ni,Ni);
Jex=[];%%%Jxm(1)*binornd(1,Px(1),Ne,Nx);
Jix=[];%%%Jxm(2)*binornd(1,Px(2),Ni,Nx);
tGen=toc;
disp(sprintf('\nTime to generate connections: %.2f sec',tGen))

%%%
% Spike times of ffwd layer
% nspikeX=poissrnd(Nx*rx*T);
% st=rand(nspikeX,1)*T;
% sx=zeros(2,numel(st));
% sx(1,:)=sort(st);
% sx(2,:)=randi(Nx,1,numel(st)); % neuron indices
% clear st;
sx=[];
nspikeX=0;
%%%

tGenx=toc;
disp(sprintf('\nTime to generate ffwd spikes: %.2f sec',tGenx))

% Neuron parameters
taum=15;
EL=-72;
Vth=0;
Vre=-75;
DeltaT=2;
VT=-55;

% Plasticity params
etaJe=3; % Learning rates
etaJi=2;
tauSTDP=200;
r0e=0.008; % Target rates
r0i=0.010;

% Random initial voltages
V0=rand(N,1)*(VT-Vre)+Vre;

% Maximum number of spikes for all neurons
% in simulation. Make it 50Hz across all neurons
% If there are more spikes, the simulation will
% terminate
maxns=ceil(.05*Ne*T);

% Time discretization of recordings. This can be coarser 
% than the dt for the Euler solver. If so, it will record
% the average current across each coarser time bin
dtRecord=500;
nBinsRecord=round(dtRecord/dt);
timeRecord=dtRecord:dtRecord:T;
Ntrec=numel(timeRecord);

% Indices of E and I neurons from which to record currents, voltages
Ierecord=1:Ne;
numerecord=numel(Ierecord);
Iirecord=1:Ni;
numirecord=numel(Iirecord);

% Integer division function
IntDivide=@(n,k)(floor((n-1)/k)+1);

Ve=V0(1:Ne);
Vi=V0(Ne+1:N);
Iee=zeros(Ne,1);
Iei=zeros(Ne,1);
Iie=zeros(Ni,1);
Iii=zeros(Ni,1);
Iex=zeros(Ne,1);
Iix=zeros(Ni,1);
xe=zeros(Ne,1);
xi=zeros(Ni,1);
VeRec=zeros(numerecord,Ntrec);
ViRec=zeros(numirecord,Ntrec);
IeeRec=zeros(numerecord,Ntrec);
IeiRec=zeros(numerecord,Ntrec);
IieRec=zeros(numirecord,Ntrec);
IiiRec=zeros(numirecord,Ntrec);
IexRec=zeros(numerecord,Ntrec);
IixRec=zeros(numirecord,Ntrec);
iXspike=1;
se=zeros(2,maxns);
si=zeros(2,maxns);
nespike=0;
nispike=0;
TooManySpikes=0;
tic
for i=1:numel(time)

    % Store recorded variables
    ii=IntDivide(i,nBinsRecord); 
    IeeRec(:,ii)=IeeRec(:,ii)+Iee(Ierecord);
    IeiRec(:,ii)=IeiRec(:,ii)+Iei(Ierecord);
    IieRec(:,ii)=IieRec(:,ii)+Iie(Iirecord);
    IiiRec(:,ii)=IiiRec(:,ii)+Iii(Iirecord);
    IexRec(:,ii)=IexRec(:,ii)+Iex(Ierecord)+Xe(Ierecord);
    IixRec(:,ii)=IixRec(:,ii)+Iix(Iirecord)+Xi(Iirecord);    
    VeRec(:,ii)=VeRec(:,ii)+Ve(Ierecord);
    ViRec(:,ii)=ViRec(:,ii)+Vi(Iirecord);
    
    
    % Propogate ffwd spikes
    %%% added: ~isempty(sx) && 
    while(~isempty(sx) && sx(1,iXspike)<=time(i) && iXspike<nspikeX)
        jpre=sx(2,iXspike);
        Iex=Iex+Jex(:,jpre)/taux;
        Iix=Iix+Jix(:,jpre)/taux;
        iXspike=iXspike+1;
    end
    
    
    % Euler update to V
    Ve=Ve+(dt/taum)*(Iee+Iei+Iex+Xe+(EL-Ve)+DeltaT*exp((Ve-VT)/DeltaT));
    Vi=Vi+(dt/taum)*(Iie+Iii+Iix+Xi+(EL-Vi)+DeltaT*exp((Vi-VT)/DeltaT));
    
    % Find which neurons spiked
     
    % If there are e spikes
    Ispike=find(Ve>=Vth);    
    if(~isempty(Ispike))

        % Store spike times and neuron indices
        if(nespike+numel(Ispike)<=maxns)
            se(1,nespike+1:nespike+numel(Ispike))=time(i);
            se(2,nespike+1:nespike+numel(Ispike))=Ispike;
        else
            TooManySpikes=1;
            break;
        end

        % Reset e mem pot.
        Ve(Ispike)=Vre;        
        
        % Update exc synaptic currents
        Iee=Iee+sum(Jee(:,Ispike),2)/taue;
        Iie=Iie+sum(Jie(:,Ispike),2)/taue;
        
        % Update cumulative number of e spikes
        nespike=nespike+numel(Ispike);
        
        % If there is plasticity onto e neurons
        if(etaJe~=0)
            Jei(Ispike,:)=Jei(Ispike,:)+Jei(Ispike,:).*repmat(etaJe*xi',numel(Ispike),1);
        end
        
        % Update rate estimates for plasticity rules
        xe(Ispike)=xe(Ispike)+1/tauSTDP;
 
    end    
    

    % If there are i spikes
    Ispike=find(Vi>=Vth);
    if(~isempty(Ispike))

        % Store spike times and neuron indices
        if(nispike+numel(Ispike)<=maxns)
            si(1,nispike+1:nispike+numel(Ispike))=time(i);
            si(2,nispike+1:nispike+numel(Ispike))=Ispike;
        else
            TooManySpikes=1;
            break;
        end

        % Reset i mem pot.
        Vi(Ispike)=Vre;        
        
        % Update inh synaptic currents
        Iei=Iei+sum(Jei(:,Ispike),2)/taui;
        Iii=Iii+sum(Jii(:,Ispike),2)/taui;
        
        % Update cumulative number of i spikes
        nispike=nispike+numel(Ispike);
        
        % If there is plasticity onto i neurons
        if(etaJi~=0)
            Jii(Ispike,:)=Jii(Ispike,:)+Jii(Ispike,:).*repmat(etaJi*xi',numel(Ispike),1);            
            Jii(:,Ispike)=Jii(:,Ispike)+Jii(:,Ispike).*repmat(etaJi*(xi-2*r0i),1,numel(Ispike));
        end
        if(etaJe~=0)
            Jei(:,Ispike)=Jei(:,Ispike)+Jei(:,Ispike).*repmat(etaJe*(xe-2*r0e),1,numel(Ispike));
        end
                    
        
        % Update rate estimates for plasticity rules
        xi(Ispike)=xi(Ispike)+1/tauSTDP;
       
    end
    
    % Euler update to synaptic currents
    Iee=Iee-dt*Iee/taue;
    Iei=Iei-dt*Iei/taui;
    Iex=Iex-dt*Iex/taux;
    Iie=Iie-dt*Iie/taue;
    Iii=Iii-dt*Iii/taui;
    Iix=Iix-dt*Iix/taux;
    
    % Update time-dependent firing rates for plasticity
    xe=xe-dt*xe/tauSTDP;
    xi=xi-dt*xi/tauSTDP;       
    
end
IeeRec=IeeRec/nBinsRecord;
IeiRec=IeiRec/nBinsRecord;
IieRec=IieRec/nBinsRecord;
IiiRec=IiiRec/nBinsRecord;
IexRec=IexRec/nBinsRecord;
IixRec=IixRec/nBinsRecord;   
VeRec=VeRec/nBinsRecord;
ViRec=ViRec/nBinsRecord;
se=se(:,1:nespike); % Get rid of padding in s
si=si(:,1:nispike); % Get rid of padding in s
tSim=toc;
disp(sprintf('\nTime for simulation: %.2f min',tSim/60))

% Compute time-dependent rates
dtRate=100;
RateTime=(dtRate/2):dtRate:T-dtRate/2;
ret=1000*hist(se(1,:),RateTime)/(dtRate*Ne);
rit=1000*hist(si(1,:),RateTime)/(dtRate*Ni);

figure
plot(RateTime,ret,'r','linewidth',2)
hold on
plot(RateTime,rit,'b','linewidth',2)
title('Time-dependent rates of ISP on EIF model')
xlabel('time (ms)')
ylabel('rate (Hz)')
legend('exc. rate', 'inh. rate')
set(gca,'FontSize',16)








