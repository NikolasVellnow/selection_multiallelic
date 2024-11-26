% DetMultiRandA.m

% Program is for n alleles.  
% It generates or uses a single random A matrix. 
% It then produces many deterministic allele frequency trajectories
% with random initial starting frequencies.
% Trajectories are coloured according to their final destination.

clear vars
close all
clc

% Parameters
n=3;                                        % No. alleles
T=1e3;                                      % Final time
s0=1e-1;                                    % Measure of absolute strength selection
Reps=1e5;

% Useful
Tset=(0:T)';
F=ones(n,1);

% Random choice of A
A=randn(n);
A=s0*(A+A');

% Initialisation
X=zeros(n,T+1,Reps);                           

% Determine trajectories
for r=1:Reps
    X0=rand(n,1); X0=X0/sum(X0);                % Random initial frequencies
    X(:,1,r)=X0;
    % time loop
    for k=1:T
        x=X(:,k,r);
        V=diag(x)-x*x';
        xp=x+V*A*x/(1+x'*A*x);
        xp=xp/sum(xp);                          % forces normalisation (may not be necessary)
        X(:,k+1,r)=xp;
    end
end

X=round(1e5*X)/1e5;                             % rounds X to 1e-5 accuracy (makes final state definite) 
a=squeeze(X(:,end,:));
u=unique(a','rows')';                           % columns are the unique end frequencies
Nsolns=size(u,2);                               % Number of unique endpoints

% Index that labels endpoint
I=zeros(Reps,1);
for r=1:Reps
    x=squeeze(X(:,:,r));
    xend=x(:,end);
    v=0;
    for i=1:Nsolns
        v=v+i*prod(xend==u(:,i));
    end
    I(r)=v;
end

c=['b','r','g'];                % colours of three possible endpoints
col=c(I);

lw=2/10;                            % linewidth trajectories
for r=1:Reps
    x=squeeze(X(:,:,r));
    plot3(x(1,:),x(2,:),x(3,:),col(r),'linewidth',lw)
    hold on
end

box("on")
set(gca,'BoxStyle','full')
grid
axis equal
axis([0,1,0,1,0,1])

% Edges around simplex
lw=3;                                   % linewidth edges
x=linspace(0, 1, 2000);
y=1-x;
plot3(x,y,0*x,'k','linewidth',lw)
plot3(x,0*x,y,'k','linewidth',lw)
plot3(0*x,x,y,'k','linewidth',lw)
plot3(x,0*x,0*x,'k','linewidth',lw)
plot3(0*x,x,0*x,'k','linewidth',lw)
plot3(0*x,0*x,x,'k','linewidth',lw)

% Axes formatting
fs=30;
set(gca,'linewidth',1,'fontsize',fs-15)
set(gca,'xtick',0:0.333:1,'xticklabel',{'0','1/3','2/3','1'})
set(gca,'ytick',0:0.333:1,'yticklabel',{'0','1/3','2/3','1'})
set(gca,'ztick',0:0.333:1,'zticklabel',{'0','1/3','2/3','1'})
set(gca,"CameraPosition",[7.3672e+00 -3.9442e+00 3.3445e+00])
xlabel('$x_1$','fontsize',fs,'Interpreter','latex');
ylabel('$x_2$','fontsize',fs,'Interpreter','latex');
zlabel('$x_3$','fontsize',fs,'Interpreter','latex');
ounits = get(gcf,'Units');
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1],'Units',ounits)
orient landscape

print -dpng n_is_3