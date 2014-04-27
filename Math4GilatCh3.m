%% CHAPTER 3 PROBLEMS: Mathmatical Operations with Arrays
%
% Brett Davis
% Chapter 3 Problems
%
% numbers: 1,4,5,6,9,10,11,14,15,17,18

%% Problem 1
clear all
close all
% evaluate the following function y for the 
% following values:

x = [-2:5]

y = (((2.*x.^2)-(5.*x)+4).^3)./x.^2

%% Problem 4
clear all
close all
% If a basketball is dropped down from a helicopter;
% its velocity as a function of time can be modeled
% by the following equation:

%vars

g = 9.81 %m/s^2;

C_d = 0.5 %kg/m^3

p = 1.2 %kg/m^3 

m = 0.624 %kg

r = 0.117 %m

ballArea = pi *( r^2)

t = [0:10] %seconds

%formula

v_of_t = sqrt((2*m*g)/(p*ballArea*C_d)) * (1 - exp((-sqrt(p*g*C_d*ballArea)/(2*m)).*t))

%% Problem 5
clear all
close all
% The magnitude of a vector u = [xi+yj+zk] is given by
% abs(u) = 14i+25j-10k ; determine its length in two
% ways:

% A) define the vector in MATLAB, then write an expression
% that uses the components of the vector:

u = [14 25 -10]

len_u = sqrt((u(1)^2) + ((u(2)^2)) + ((u(3)^2)))

% B) define the vector in matlab then use element by 
% element operation to create a new vector with elements
% that are the square of the original vector.

U = u.^2

len_U = sqrt(sum(U))

%% Problet 6
clear all
close all
% The position as a function of time (x(t),y(t)) of a 
% projectile fired with a speed of v_0 at an angle theta
% is given by some formulas, given some data use the 
% provided formulas to calcuate the distance r for:

% vars

t = [0:2:20] %secs

g = 9.81 %m/s^2

v_0 = 100 %m/s

theta = 79 %degs

% formulas

x_of_t = (v_0 * cosd(theta) .*t)

y_of_t = (v_0 * sind(theta) .*t) - ((1/2)*g.*t.^2)
 
r_of_t = sqrt((x_of_t.^2)+(y_of_t.^2))

%% Problem 9
clear all
close all
% Define h and k as the following scalars:

h = 0.7
k = 8.85

% and x y and z as the following vectors:

x = [1:5]
y = [2.1:-.1:1.7]
z = [2.0:0.5:4.0]

% use these elements to calculate G using element-by
% element calculations for the vectors.

G = (((h.*x)+(k.*y))./(x+y).^h) + (exp((h.*y)./z))./(z.^(y./x))

%% Problem 10
clear all
close all
format long
% show that the limit as x approaches zero of the 
% following function is equal to 1

%vars

x = [1 .5 .1 .01 .001 .00001 .0000001]

y = (exp(x)-1)./x

%% Problem 11
clear all
close all
format short
% Use matlab to show that the sum of the infinite series
% given converes to pi

%vars

n_100 = [0:100]

n_10k = [0:10000];

n_mil = [0:1000000];

%formula

a = sum(4*((-1).^n_100)./((2.*n_100)+1))

b = sum(4*((-1).^n_10k)./((2.*n_10k)+1))

c = sum(4*((-1).^n_mil)./((2.*n_mil)+1))

%% Problem 14
clear all
close all
format short
% Make the following matrices

A = [5 2 4;1 7 -3;6 -10 0]

B = [11 5 -3; 0 -12 4; 2 6 1]

C = [7 14 1;10 3 -2;8 -5 9]

% A) calculate A + B and B + A to show that addition 
% of matrices is commutative

sum_1 = A + B

sum_2 = B + A

% B) Cacluate A + (B+C) and (A+B) + C

sum_3 = A + (B+C)

sum_4 = (A+B) + C

% C) Calculate:

sum_5 = 5*(A+C)

sum_6 = (5*A) + (5*C)

% D) Calculate: 

sum_7 = A * (B + C)

sum_8 = (A*B) + (A*C)

%% Problem 15

% Using the matrices from the previous problem to answer the
% following:

% A) Does

a_1 = A * B

% equal:

a_2 = B * A

% No, it does not!

% B) Does

b_1 = A*(B*C)

% equal:

b_2 = (A*B)*C

% Yes, it does!

% C) Does

c_1 = (A*B)'

% equal:

c_2  = B' * A'

% Yes, it does!

% D) Does

d_1 = (A+B)'

% equal:

d_2 = A' + B'

% Yes, it does!

%% Problem 17
clear all
close all
format long
% The mechanical power output in a contracting muscle is
% given by a formula. T is muscle tenstion, v is the 
% shortening velocity (max of V_max) T_0 is the isometric
% tension (The tension at zero velocity) and k is a non-
% dimensional constant that ranges between .15 and .25 for
% most muscles.

% A) create the following vector u

u = [0:.05:1]

% k is equal to .25

k = .25

% B) calculate the value of P for each value of u using the 
% following formula.


p = ((k.*u) .* (1-u)) ./ (k+u)

% C) using the max function find the maximum value of p

p_max = max(p)

% D) repeat the first three steps with increments of .01 and
% find the percent relative error with a given formula.

u_small = [0:.01:1]

p_small = ((k.*u_small) .* (1-u_small)) ./ (k+u_small)

p_small_max = max(p_small)

% formula

relativeError = abs((p_small_max - p_max)/p_max) * 100

%% Problem 18
clear all
close all
format short
% solve the following system of equations:

trix_A = [1.5 -2 1 3 .5 7.5;3 1 -1 4 -3 16; 2 6 -3 -1 3 78;5 2 4 -2 6 71;-3 3 2 5 4 54]

trix_A_solved = rref(trix_A)


