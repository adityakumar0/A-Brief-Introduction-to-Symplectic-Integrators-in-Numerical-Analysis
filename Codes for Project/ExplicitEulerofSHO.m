%Code for Explicit Euler of Simple Harmonic Oscillator(S.H.O)
%t_exact is the time scale for exact solution to S.H.O
%y_exact is the exact solution of S.H.O

clear

y = ForwardEuler('harmosc',[1;0],300,0.1);
t_exact = 0:0.1:30;
y_exact = [cos(t_exact);-sin(t_exact)]; %Analtic solution of S.H.O
Eex=(y_exact(1,:)).^2+(y_exact(2,:)).^2; %Exact Energy of S.H.O
Enu=(y(1,:)).^2+(y(2,:)).^2; %Energy Computed via Explicit Euler of S.H.O

disp(Eex)
disp(Enu)

%Error
abserror=abs(y_exact - y);
disp(abserror)

%Plots
figure(1);
plot(t_exact,y_exact,0:0.1:30, y,'-o');
title('Plot of Analytic vs Numerical Solutions',... 
  'FontWeight','bold')
xlabel('Time')
ylabel('Amplitude')

figure(2);
plot(y_exact(1,:), y_exact(2,:), y(1,:), y(2,:), '-o')
title('The Phase Plane Portrait',... 
  'FontWeight','bold')

figure(3);
plot(t_exact,Eex,t_exact,Enu,'-o')
title('Plot of Analytic vs Numerical Energy',... 
  'FontWeight','bold')
xlabel('Time')
ylabel('2*Energy')