% DetMultiRandA.m

% Program is for n=3 alleles.  
% It uses a single A matrix. 
% It then produces Reps sets of random initial starting frequencies
% and from these produces deterministic allele frequency trajectories.
% Final point is black.

clear vars
close all
clc

% Parameters
n=3;                                        % No. alleles
T=1e4;                                      % Final time (in unit time steps)
Reps=20000;
A =[0.0108    0.0270   -0.0269
    0.0270    0.0064   -0.0097
   -0.0269   -0.0097    0.0716];

% Useful
I=eye(n);
Tset=(0:T)';
F=ones(n,1);

% Initialisation
X=zeros(n,T+1,Reps);                   % X(i,t,r)=frequency of allele i at time t-1 in replicate r
for r=1:Reps
    X0=rand(n,1); X0=X0/sum(X0);        
    X(:,1,r)=X0;                        % Assign random initial frequencies
    % Time loop
    for k=1:T
        x=X(:,k,r);
        V=diag(x)-x*x';
        xp=x+V*A*x;
        xp=xp/sum(xp);                  % Ensuring normalisation
        X(:,k+1,r)=xp;
    end
end

p=squeeze(abs(X(3,end,:)-1)<1e-2);      % Indicates trajectories where X(3) goes to 1
q=squeeze(abs(X(2,end,:)-0.44)<1e-2);   % Indicates trajectories where X(2) goes to 0.44

X1=X(:,:,p);                            % Set of trajectories where X(3) goes to 1               
X2=X(:,:,q);                            % Set of trajectories where X(2) goes to 0.44      

hold on
for r=1:sum(p)
    p1=plot3(X1(1,:,r),X1(2,:,r),X1(3,:,r),'b','linewidth',1/10);
    p1.Color(4)=0.15;
end
plot3(X1(1,end,r),X1(2,end,r),X1(3,end,r),'k.','MarkerSize',40)
for r=1:sum(q)
    p1=plot3(X2(1,:,r),X2(2,:,r),X2(3,:,r),'r','linewidth',1/10);
    p1.Color(4)=0.15;
end
d=0.01;
plot3(X2(1,end,r)+d,X2(2,end,r)+d,X2(3,end,r)+d,'k.','MarkerSize',40)

% B=squeeze(X1(:,1,:));
% plot3(X1(1,:),X1(2,:),X1(3,:),'b.')

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
print -dpng SpecialA