%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
% Algorithm 2: First NSFD Method For Ethanol Metabolism %
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

X2      = zeros(3,length(t));
eig2      = zeros(3,length(t));
A_0     = 1.00;
B_0     = 0.10;
C_0     = 0.01;
eig2(1,1) = 1.00;
eig2(2,1) = 1.00;
eig2(3,1) = 1.00;
X2(1,1) = A_0;
X2(2,1) = B_0;
X2(3,1) = C_0;

% Step 2.3: Eigenvalues and Eigenvectors of System Matrix of Super-Solution

A_Super   = [-a 0 0; a -b b; 0 b -(b+c/(d+A_0+B_0+C_0))];
[V,W]     = eig(A_Super);

% Step 3: NSFD Method For Solving The Differential Equation

% Step 3.1: Definition Of Matrix A2

A2      = zeros(3,3);
A2(1,1) = 1/(1+h*a);

% Step 3.2: Solution Loop For NSFD Scheme

for j = 1:1:(length(t)-1)
  A2(2,1)     = ((h*a)/((1+h*a)*(1+h*b)));
  A2(2,2)     = (1/(1+h*b));
  A2(2,3)     = (h*b)/((1+h*b));
  A2(3,1)     = (h^2*a*b)/((1+h*a)*(1+h*b)*(1+(h*c)/(d+X2(3,j))));
  A2(3,2)     = (h*b)/((1+h*b)*(1+(h*c)/(d+X2(3,j))));
  A2(3,3)     = 1/((1+h*b)*(1+(h*c)/(d+X2(3,j))));
  eig2(:,j+1) = eig(A2);
  X2(:,j+1)   = A2*X2(:,j);
endfor

% Step 4: Plotting Of Solutions

A2_max = max(max(A2(:,2:1:end)));

figure(2)

hold on

subplot(2,2,1);
plot(t,X2(1,:),'color','blue','linewidth',1.25);
hold on
plot(t,(max([A_0,B_0,C_0]))*exp(max(min(W))*t),'color','red','linewidth',1.25);
title('Amount of ethanol in stomach','fontsize',14);
xlabel('t','fontsize',12);
ylabel('A(t)','fontsize',12);
yticks([0 0.2 0.4 0.6 0.8 1.0 1.2]);

subplot(2,2,2);
plot(t,X2(2,:),'color','blue','linewidth',1.25);
hold on
plot(t,(max([A_0,B_0,C_0]))*exp(max(min(W))*t),'color','red','linewidth',1.25);
title('Amount of ethanol in blood','fontsize',14);
xlabel('t','fontsize',12);
ylabel('B(t)','fontsize',12);
yticks([0 0.2 0.4 0.6 0.8 1.0 1.2]);

subplot(2,2,3);
plot(t,X2(3,:),'color','blue','linewidth',1.25);
hold on
plot(t,(max([A_0,B_0,C_0]))*exp(max(min(W))*t),'color','red','linewidth',1.25);
title('Amount of ethanol in liver','fontsize',14);
xlabel('t','fontsize',12);
ylabel('C(t)','fontsize',12);
yticks([0 0.2 0.4 0.6 0.8 1.0 1.2]);

subplot(2,2,4);
plot(t,X2(1,:)+X2(2,:)+X2(3,:),'color','blue','linewidth',1.25);
hold on
plot(t,3*(max([A_0,B_0,C_0]))*exp(max(min(W))*t),'color','red','linewidth',1.25);
title('Total amount of ethanol in body','fontsize',14);
xlabel('t','fontsize',12);
ylabel('N(t)','fontsize',12);
yticks([0 1 2 3]);

hold off
