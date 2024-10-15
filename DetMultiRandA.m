% DetMultiRandA.m

% Program is for n=3 alleles.  
% It generates a single random A matrix. 
% It then produces 50 deterministic allele frequency trajectories
% that start from random initial starting frequencies.
% Starting points are blue, final point is black.

% Run the program a few times.
% Sometimes all trajectories go to same point - which may be equiv to
% fixation or perhaps a polymorphism, sometimes some trajectories go to one
% point, while others go to another. 

clear vars
close all
clc

% Parameters
n=3;                                        % No. alleles
T=1e4;                                      % Final time
dt=1;                                       % Time step
s0=1e-2;                                    % Measure of absolute strength selection
Reps=50;

% Calculated
Nt=T/dt;                                    % No. time steps

% Useful
I=eye(n);
Tset=(0:Nt)'*dt;
F=ones(n,1);

% Random choice of A
A=randn(n);
A=s0*(A+A');

% Initialisation
X=zeros(n,Nt+1,Reps);                           % Random starting frequencies

for r=1:Reps
    X0=rand(n,1); X0=X0/sum(X0);                % Initial frequencies
    X(:,1,r)=X0;
    % Loop
    for k=1:Nt
        x=X(:,k,r);
        V=diag(x)-x*x';
        xp=x+dt*V*A*x;
        xp=xp/sum(xp);
        X(:,k+1,r)=xp;
    end
end

for r=1:Reps
    plot3(X(1,1,r),X(2,1,r),X(3,1,r),'b.','MarkerSize',20)
    hold on
    plot3(X(1,end,r),X(2,end,r),X(3,end,r),'k.','MarkerSize',30)
    x=squeeze(X(:,:,r));
    plot3(x(1,:),x(2,:),x(3,:),'m','linewidth',1)
end

box("on")
set(gca,'BoxStyle','full')
grid
axis equal
axis([0,1,0,1,0,1])

fs=30;
set(gca,'linewidth',3,'fontsize',fs-15)
set(gca,'xtick',0:0.333:1,'xticklabel',{'0','1/3','2/3','1'})
set(gca,'ytick',0:0.333:1,'yticklabel',{'0','1/3','2/3','1'})
set(gca,'ztick',0:0.333:1,'zticklabel',{'0','1/3','2/3','1'})
set(gca,"CameraPosition",[8.4352 -2.2037 2.6731])
xlabel('$x_1$','fontsize',fs,'Interpreter','latex');
ylabel('$x_2$','fontsize',fs,'Interpreter','latex');
zlabel('$x_3$','fontsize',fs,'Interpreter','latex');
ounits = get(gcf,'Units');
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1],'Units',ounits)
orient landscape
%print -dpng SingleTrajHHP