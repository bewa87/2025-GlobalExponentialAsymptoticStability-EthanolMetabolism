%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
% Algorithm 1: First NSFD Method For Ethanol Metabolism %
%                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 0: Clearing

clear all;
close all;
clc;

% Step 1: Define All Constant Problem Parameters

a = 0.1;
b = 0.025;
c = 0.020;
d = 1;

% Step 2: Set Up Time Vector, Solution Vectors And Initial Conditions

% Step 2.1: Definition Of Time Vector

T = 1000;
h = 18;
t = 0:h:T;

% Step 2.2: Definition Of Initial Conditions And Solution Matrix X

X1        = zeros(3,length(t));
eig1      = zeros(3,length(t));
A_0       = 1.00;
B_0       = 0.10;
C_0       = 0.01;
eig1(1,1) = 1.00;
eig1(2,1) = 1.00;
eig1(3,1) = 1.00;
X1(1,1)   = A_0;
X1(2,1)   = B_0;
X1(3,1)   = C_0;

% Step 2.3: Eigenvalues and Eigenvectors of System Matrix of Super-Solution

A_Super   = [-a 0 0; a -b b; 0 b -(b+c/(d+A_0+B_0+C_0))];
[V,W]     = eig(A_Super);

% Step 3: NSFD Method For Solving The Differential Equation

% Step 3.1: Definition Of Matrix A1

A1      = zeros(3,3);
A1(1,1) = 1/(1+h*a);

% Step 3.2: Solution Loop For NSFD Scheme

for j = 1:1:(length(t)-1)
  A1(2,1)     = ((h*a)/((1+h*a)*(1+h*b)))...
                + (h^3*a*b^2)/((1+h*a)*(1+h*b)^2*(1+(h*b)/(1+h*b)+(h*c)/(d+X1(3,j))));
  A1(2,2)     = (1/(1+h*b))...
                + (h^2*b^2)/((1+h*b)^2*(1+(h*b)/(1+h*b)+(h*c)/(d+X1(3,j))));
  A1(2,3)     = (h*b)/((1+h*b)*(1+(h*b)/(1+h*b)+(h*c)/(d+X1(3,j))));
  A1(3,1)     = (h^2*a*b)/((1+h*a)*(1+h*b)*(1+(h*b)/(1+h*b)+(h*c)/(d+X1(3,j))));
  A1(3,2)     = (h*b)/((1+h*b)*(1+(h*b)/(1+h*b)+(h*c)/(d+X1(3,j))));
  A1(3,3)     = 1/(1+(h*b)/(1+h*b)+(h*c)/(d+X1(3,j)));
  eig1(:,j+1) = eig(A1);
  X1(:,j+1)   = A1*X1(:,j);
endfor

% Step 4: Plotting Of Solutions

A1_max = max(max(A1(:,2:1:end)));

figure(1)

hold on

subplot(2,2,1);
plot(t,X1(1,:),'color','blue','linewidth',1.25);
hold on
plot(t,(max([A_0,B_0,C_0]))*exp(max(min(W))*t),'color','red','linewidth',1.25);
title('Amount of ethanol in stomach','fontsize',14);
xlabel('t','fontsize',12);
ylabel('A(t)','fontsize',12);
yticks([0 0.2 0.4 0.6 0.8 1.0 1.2]);
hold off

subplot(2,2,2);
plot(t,X1(2,:),'color','blue','linewidth',1.25);
hold on
plot(t,(max([A_0,B_0,C_0]))*exp(max(min(W))*t),'color','red','linewidth',1.25);
title('Amount of ethanol in blood','fontsize',14);
xlabel('t','fontsize',12);
ylabel('B(t)','fontsize',12);
yticks([0 0.2 0.4 0.6 0.8 1.0 1.2]);

subplot(2,2,3);
plot(t,X1(3,:),'color','blue','linewidth',1.25);
hold on
plot(t,(max([A_0,B_0,C_0]))*exp(max(min(W))*t),'color','red','linewidth',1.25);
title('Amount of ethanol in liver','fontsize',14);
xlabel('t','fontsize',12);
ylabel('C(t)','fontsize',12);
yticks([0 0.2 0.4 0.6 0.8 1.0 1.2]);

subplot(2,2,4);
plot(t,X1(1,:)+X1(2,:)+X1(3,:),'color','blue','linewidth',1.25);
hold on
plot(t,3*(max([A_0,B_0,C_0]))*exp(max(min(W))*t),'color','red','linewidth',1.25);
title('Total amount of ethanol in body','fontsize',14);
xlabel('t','fontsize',12);
ylabel('N(t)','fontsize',12);
yticks([0 1 2 3]);

hold off
