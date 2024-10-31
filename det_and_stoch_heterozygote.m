% Case of calculation of fitness effects with heterozygote advantage

% Program is for n=3 alleles.  
% It generates a single 'A' matrix for the case with heterozygote advantage. 
% It produces a total of 'Reps' deterministic allele frequency trajectories
% that start from starting frequencies that have to be provided.
% Starting points are blue, final points (of deterministic trajectories) are red. 
% Deterministic trajectories are in magenta.
% It also produces a total of 'Reps' stochastic allele frequency trajectories
% that start from same initial starting frequencies as deterministic trajectories.
% Stochastic trajectories are in black.

clear vars
close all
clc

% Parameters
n=3;                                        % No. alleles
T=2e4;                                      % Final time
sigma = 0.05;                               % advantage of heterozygotes
N=5e4;                                      % Population size
Reps=3;                                     % No. deterministic trajectories

% Create matrix of fitness effects A
A = calc_A_heterozygote(sigma, n);

% Initialisation
Xd=zeros(n,T+1,Reps);                           % Stores deterministic trajectories
Xs=zeros(n,T+1,Reps);                           % Stores stochastic trajectories


for r=1:Reps
    % Choose starting frequencies for first three replicates
    if r==1
        X0=[0.7; 0.15; 0.15]; X0=X0/sum(X0);    % Forces normalisation (just to be sure!)
    elseif r==2
        X0=[0.15; 0.7; 0.15]; X0=X0/sum(X0);    % Forces normalisation (just to be sure!)
    elseif r==3
        X0=[0.15; 0.15; 0.7]; X0=X0/sum(X0);    % Forces normalisation (just to be sure!)
    else
        X0=rand(n,1); X0=X0/sum(X0);
    end
    Xd(:,1,r)=X0;
    Xs(:,1,r)=X0;
    % Deterministic loop
    for k=1:T
        x=Xd(:,k,r);
        V=diag(x)-x*x';
        D=V*A*x/(1+x'*A*x);
        xp=x+D;
        xp=xp/sum(xp);                          % Forces normalisation (just to be sure!)
        Xd(:,k+1,r)=xp;
    end
    % Stochastic loop
    for k=1:T
        x=Xs(:,k,r);
        V=diag(x)-x*x';
        D=V*A*x/(1+x'*A*x);
        m=mnrnd(2*N,x+D)';                      % Multinomial random variable (implements drift)
        xp=m/2/N;
        Xs(:,k+1,r)=xp;
    end
end

% plotting

for r=1:Reps
    plot3(Xd(1,1,r),Xd(2,1,r),Xd(3,1,r),'b.','MarkerSize',25)               % plots initial point of all trajectories
    hold on
    plot3(Xd(1,end,r),Xd(2,end,r),Xd(3,end,r),'r.','MarkerSize',40)         % plots final point(s) deterministic trajectories
    x=squeeze(Xd(:,:,r));
    plot3(x(1,:),x(2,:),x(3,:),'m','linewidth',1)                           % plots deterministic trajectories
    x=squeeze(Xs(:,:,r));
    plot3(x(1,1:T/10),x(2,1:T/10),x(3,1:T/10),'k','linewidth',1)                           % plots stochastic trajectories
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

%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the grid for x and y
[x, y] = meshgrid(linspace(0, 1, 2000), linspace(0, 1, 2000));

% Calculate z based on the plane equation x + y + z = 1
z=1-x-y;

% Set z values to NaN where they are less than 0 to avoid plotting outside the unit cube
z(z < 0) = NaN;

% Plot the surface of the plane with a solid color (e.g., blue) and partial transparency
surf(x, y, z, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.1);

% % Set axis limits
% xlim([0 1]);
% ylim([0 1]);
% zlim([0 1]);
% 
% % Add labels and title
% xlabel('X-axis');
% ylabel('Y-axis');
% zlabel('Z-axis');
% title('Partly Transparent Shaded Plane in Unit Cube');
% 
% % Add grid for better visualization
% grid on;
% view(3); % Set a 3D view
% hold off;
% axis equal
% box on
% set(gca,'boxstyle','full')

%%%%%%%%%%%%%%%%%%%%%%%%%%

ounits = get(gcf,'Units');
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1],'Units',ounits)
orient landscape
print -dpng Both



