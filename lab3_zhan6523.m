%% ODE Lab: Creating your own ODE solver in MATLAB
%
% In this lab, you will write your own ODE solver for the Improved Euler 
% method (also known as the Heun method), and compare its results to those 
% of |ode45|.
%
% You will also learn how to write a function in a separate m-file and 
% execute it.
% 
% Opening the m-file lab3.m in the MATLAB editor, step through each
% part using cell mode to see the results.  Compare the output with the
% PDF, which was generated from this m-file.
%
% There are six (6) exercises in this lab that are to be handed in on the
% due date. Write your solutions in the template, including
% appropriate descriptions in each step. Save the .m files and submit them 
% online on Quercus.
%
% MAT292, Fall 2019, Stinchcombe & Parsch, modified from
% MAT292, Fall 2018, Stinchcombe & Khovanskii, modified from
% MAT292, Fall 2017, Stinchcombe & Sinnamon, modified from
% MAT292, Fall 2015, Sousa, modified from
% MAT292, Fall 2013, Sinnamon & Sousa, modified from
% MAT292, Fall 2011, Hart & Pym

%% Student Information
%
% Student Name: Xiao En (Shawn) Zhang
%
% Student Number: 1005154921
%

%% Creating new functions using m-files.
%  
% Create a new function in a separate m-file:
%
% Specifics:  Create a text file with the file name f.m
% with the following lines of code (text):
%
%  function y = f(a,b,c) 
%  y = a+b+c;
%
% Now MATLAB can call the new function f (which simply accepts 3 numbers
% and adds them together).  
% To see how this works, type the following in the matlab command window:
% sum = f(1,2,3)

%% Exercise 1
%
% Objective: Write your own ODE solver (using the Heun/Improved Euler
% Method).
%
% Details: This m-file should be a function which accepts as variables 
% (t0,tN,y0,h), where t0 and tN are the start and end points of the 
% interval on which to solve the ODE, y0 is the initial condition of the
% ODE, and h is the stepsize.  You may also want to pass the function into
% the ODE the way |ode45| does (check lab 2).
%
% Note: you will need to use a loop to do this exercise.  
% You will also need to recall the Heun/Improved Euler algorithm learned in lectures.  

%% Exercise 2
%
% Objective: Compare Heun with |ode45|.
%
% Specifics:  For the following initial-value problems (from lab 2, 
% exercises 1, 4-6), approximate the solutions with your function from
% exercise 1 (Improved Euler Method).
% Plot the graphs of your Improved Euler Approximation with the |ode45| 
% approximation.

% (a) |y' = y tan t + sin t, y(0) = -1/2| from |t = 0| to |t = pi|
f = @(t,y) y.*tan(t) + sin(t); 

% The initial conditions
t0 = 0;
t1 = pi;
y0 = -1/2;
z1 = heun(f,t0,t1,y0,0.01);
% z=euler(f,y0,0:0.01:1)
n = (t1-t0)/0.01;
plot(linspace(t0,t1,n),z1);

soln = ode45(f, [t0, t1], y0);
plot(linspace(t0,t1,n), z1, soln.x, soln.y, 'x', 'MarkerSize',10, 'LineWidth', 2);
xlabel('t');
ylabel('y');
legend('Numerical', 'Exact','Location','Best');
% (b) |y' = 1 / y^2 , y(1) = 1| from |t=1| to |t=10|

f = @(t,y) 1./(y.^2); 

% The initial conditions
t0 = 1;
t1 = 10;
y0 = 1;

z1 = heun(f,t0,t1,y0,0.01);
% z=euler(f,y0,0:0.01:1)
n = (t1-t0)/0.01;
plot(linspace(t0,t1,n),z1);

soln = ode45(f, [t0, t1], y0);
plot(linspace(t0,t1,n), z1, soln.x, soln.y, 'x', 'MarkerSize',10, 'LineWidth', 2);
xlabel('t');
ylabel('y');
legend('Numerical', 'Exact','Location','Best');
% (c) |y' =  1 - t y / 2, y(0) = -1| from |t=0| to |t=10|
f = @(t,y) 1 - t.*y./2; 

% The initial conditions
t0 = 0;
t1 = 10;
y0 = -1;

%solution
soln = ode45(f, [t0, t1], y0);

z1 = heun(f,t0,t1,y0,0.01);
% z=euler(f,y0,0:0.01:1)
n = (t1-t0)/0.01;
plot(linspace(t0,t1,n),z1);

soln = ode45(f, [t0, t1], y0);
plot(linspace(t0,t1,n), z1, soln.x, soln.y, 'x', 'MarkerSize',10, 'LineWidth', 2);
xlabel('t');
ylabel('y');
legend('Numerical', 'Exact','Location','Best');
% (d) |y' = y^3 - t^2, y(0) = 1| from |t=0| to |t=1|
f = @(t,y) y.^3 - t.^2; 

% The initial conditions
t0 = 0;
t1 = 1;
y0 = 1;

%solution
soln = ode45(f, [t0, t1], y0);

z1 = heun(f,t0,t1,y0,0.01);
% z=euler(f,y0,0:0.01:1)
n = (t1-t0)/0.01;
plot(linspace(t0,t1,n),z1);

soln = ode45(f, [t0, t1], y0);
plot(linspace(t0,t1,n), z1, soln.x, soln.y, 'x', 'MarkerSize',10, 'LineWidth', 2);
xlabel('t');
ylabel('y');
legend('Numerical', 'Exact','Location','Best');
% Comment on any major differences, or the lack thereof. You do not need
% to reproduce all the code here. Simply make note of any differences for
% each of the four IVPs.
% For part a, there was a small discrepancy around t=0, likely due to the change in sign of the derivative.
% For part d, the numerical approximations blew up due to a vertical
% tangent in the real solution.

%% Exercise 3
%
% Objective: Use Euler's method and verify an estimate for the global error.
%
% Details: 
%
% (a) Use Euler's method (you can use
% euler.m from iode) to solve the IVP
%
% |y' = 2 t sqrt( 1 - y^2 )  ,  y(0) = 0|
f = @(t,y) 2*t*sqrt(1-y^2);

% from |t=0| to |t=0.5|.
y0 = 0;

z=euler(f,y0,0:0.02:0.5)

% (b) Calculate the solution of the IVP and evaluate it at |t=0.5|.
% Exact Solution y=sin(t^2)
% y(0.5) = 0.247403959...

% (c) Read the attached derivation of an estimate of the global error for 
%     Euler's method. Type out the resulting bound for En here in
%     a comment. Define each variable.

%En = (1 + M)(delta t)/2 * (e^(M*delta t*n) - 1)
% M = 1 is an upper bound, delta t is the time step (0.01,0.02), n is the number of
% time steps (51,26)



% (d) Compute the error estimate for |t=0.5| and compare with the actual
% error.
% En = 0.00648
% actual error = 0.0047039
% Seems fairly accurate for 50 time steps

% (e) Change the time step and compare the new error estimate with the
% actual error. Comment on how it confirms the order of Euler's method.

% new error (26 time steps): 0.01364
% actual error = 0.00950

%The order of the Euler Method is shown to be first order as the error doubles when the
%time step gets halved.
%% Adaptive Step Size
%
% As mentioned in lab 2, the step size in |ode45| is adapted to a
% specific error tolerance.
%
% The idea of adaptive step size is to change the step size |h| to a
% smaller number whenever the derivative of the solution changes quickly.
% This is done by evaluating f(t,y) and checking how it changes from one
% iteration to the next.

%% Exercise 4
%
% Objective: Create an Adaptive Euler method, with an adaptive step size |h|.
%
% Details: Create an m-file which accepts the variables |(t0,tN,y0,h)|, as 
% in exercise 1, where |h| is an initial step size. You may also want to 
% pass the function into the ODE the way |ode45| does.
%
% Create an implementation of Euler's method by modifying your solution to 
% exercise 1. Change it to include the following:
%
% (a) On each timestep, make two estimates of the value of the solution at
% the end of the timestep: |Y| from one Euler step of size |h| and |Z| 
% from two successive Euler steps of size |h/2|. The difference in these
% two values is an estimate for the error.
%
% (b) Let |tol=1e-8| and |D=Z-Y|. If |abs(D)<tol|, declare the step to be
% successful and set the new solution value to be |Z+D|. This value has
% local error |O(h^3)|. If |abs(D)>=tol|, reject this step and repeat it 
% with a new step size, from (c).
%
% (c) Update the step size as |h = 0.9*h*min(max(tol/abs(D),0.3),2)|.
%
% Comment on what the formula for updating the step size is attempting to
% achieve.
% By decreasing the step size, the local error for that time step would
% decrease.

%% Exercise 5
%
% Objective: Compare Euler to your Adaptive Euler method.
%
% Details: Consider the IVP from exercise 3.
%
% (a) Use Euler method to approximate the solution from |t=0| to |t=0.75|
% with |h=0.025|.
f = @(t,y) 2*t*sqrt(1-y^2);

% from |t=0| to |t=0.75|.
y0 = 0;
timestep = 0.01;

z=euler(f,y0,0:timestep:0.75)
a=adaptive_euler(f,0,0.75,0,timestep);
soln = ode45(f, [0, 0.75], y0);

%Plot both approximations together with the exact solution.
n=0.75/timestep + 1;
x=linspace(0,0.75,n);

plot(x, z, x, a, soln.x, soln.y, 'x', 'MarkerSize',10, 'LineWidth', 2);
xlabel('t');
ylabel('y');
legend('Euler', 'Adaptive','Exact','Location','Best');
%% Exercise 6
%
% Objective: Problems with Numerical Methods.
%
% Details: Consider the IVP from exercise 3 (and 5).
% 
% (a) From the two approximations calculated in exercise 5, which one is
% closer to the actual solution (done in 3.b)? Explain why.
% 
%The adaptive method was closer to the actual solution. This is because
%once the local error goes beyond a certain threshold, the method takes
%smaller timesteps until the error is within an acceptable error range.

% (b) Plot the exact solution (from exercise 3.b), the Euler's 
% approximation (from exercise 3.a) and the adaptive Euler's approximation 
% (from exercise 5) from |t=0| to |t=1.5|.
f = @(t,y) 2*t*sqrt(1-y^2);

y0 = 0;
t1 = 1.5;
t0 = 0;
timestep = 0.01;
n=t1/0.01 + 1;

z=euler(f,y0,0:timestep:t1)
a=adaptive_euler(f,0,t1,0,timestep);

tt = linspace(t0,t1,n);
yy = sin(tt.^2);

x=linspace(t0,t1,n);
%Plotting Graphs
plot(x, z, x, a, tt, yy);
xlabel('t');
ylabel('y');
legend('Euler', 'Adaptive','Exact','Location','Best');
% (c) Notice how the exact solution and the approximations become very
% different. Why is that? Write your answer as a comment.
% Near y=1, the slope approaches 0 but then becomes negative for the actual
% solution. For the approximations,  the timesteps estimate the next point
% using the current and estimated next point's slope, which will result in
% a positive or 0 slope. This is why the approximations form straight lines
% after the sin wave hits a maximum.
