%Code for Forward Euler
%f is the function for system of ODEs of S.H.O (x' and y')
%y0 is the initial conditions for system of ODEs
%N is the number of steps
%h is the step size

function ys = ForwardEuler(f, y0, N, h)
ys = zeros(length(y0), N+1);
y = y0;
ys(:,1) = y0;
for k = 1:N
y = y + h * feval(f, y);
ys(:,k+1) = y;
end

end


