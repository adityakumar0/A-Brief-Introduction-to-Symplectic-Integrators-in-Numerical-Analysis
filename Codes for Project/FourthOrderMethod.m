%Non-symplectic fourth order ode built-in solver used for S.H.O
%t_ex is the time scale for exact solution to S.H.O
%y_ex is the exact solution of S.H.O


clear

%Numerical Solution of S.H.O
timerange=[0:0.1:10000];
initialvalues=[1 0];
options=odeset('OutputFcn','odephas2');
ode45(@f,timerange,initialvalues,options)
[t,y]=ode45(@f,timerange,initialvalues);

t_ex=0:0.1:10000;
y_ex=[cos(t_ex);-sin(t_ex)]; %Analytic Solution of S.H.O

Eex=((y_ex(1,:)).^2+(y_ex(2,:)).^2).'; %Exact Energy of S.H.O
Enu=(y(:,1)).^2+(y(:,2)).^2; %Energy Computed via Explicit Euler of S.H.O

%Error
abserror=abs((y_ex).' - y);
disp(abserror)


%Plots
figure(1);
plot(t_ex,y_ex,t,y,'-o');
title('Plot of Analytic vs Numerical Solutions',... 
  'FontWeight','bold')
xlabel('Time')
ylabel('Amplitude')

figure(2);
plot(y_ex(1,:),y_ex(2,:),y(:,1),y(:,2))
title('The Phase Plane Portrait',... 
  'FontWeight','bold')

figure(3);
plot(t_ex,Eex,t_ex,Enu,'-o')
title('Plot of Analytic vs Numerical Energy',... 
  'FontWeight','bold')
xlabel('Time')
ylabel('2*Energy')