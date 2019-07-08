%Symplectic Explicit Euler
%p is calculated using Explicit Euler
%q is calculated using Implicit Euler
%The entire method is just a combination of the two Euler methods being
%solved simulataneously
%Symplectic Explicit Euler Method

clear

dpdt=@(q)(-q);
dqdt=@(p)(p);
p0=0;
q0=1;
tspan=[0:0.01:1000];

%Numerical Solution of S.H.O
N = length(q0);  
Nt = length(tspan); 
hs = diff(tspan);
q = zeros(N,Nt); q(:,1)=q0;
p = zeros(N,Nt); p(:,1)=p0;
for nt = 2:Nt
   h = hs(nt-1);
	p(:,nt) = p(:,nt-1) + h * feval(dpdt, q(:,nt-1)); 
   q(:,nt) = q(:,nt-1) + h * feval(dqdt, p(:,nt));
end
q = q.';
p = p.';

y_nu=[q.';p.'];
t_ex=[0:0.01:1000];
y_ex=[cos(t_ex);-sin(t_ex)]; %Analytic Solution of S.H.O

Eex=(y_ex(1,:)).^2 + (y_ex(2,:)).^2; %Exact Energy of S.H.O
Enu=(q.').^2 + (p.').^2; %Numerically Computed Energy of S.H.O

%Error
abserror=abs(y_ex - y_nu);
disp(abserror)

%Plot
figure(1);
plot(t_ex , y_ex, t_ex, y_nu);
title('Plot of Analytic vs Numerical Solutions',... 
  'FontWeight','bold')
xlabel('Time')
ylabel('Amplitude')

figure(2);
set(gca,'nextplot','replacechildren'); % axis settings
lsh=plot(y_ex(1,:),y_ex(2,:), q.', p.','Color',[0 0 1]); % plot
% axis settings
set(lsh(1),'Color',[1 0 0]);
set(lsh(2),'Color',[0 0 1]);
axis square;
axis([-1 1 -1 1]);
title('Phase space evolution; y''''+y = 0')
xlabel(gca,'q (= y)');
ylabel('p (= y'')');

figure(3);
plot(t_ex,Eex,t_ex,Enu)
title('Plot of Analytic vs Numerical Energy',... 
  'FontWeight','bold')
xlabel('Time')
ylabel('2*Energy')
