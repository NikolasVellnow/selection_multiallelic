% Case of randomly chosen genotype fitness effects

% Program is for n=3 alleles.  
% It generates a single 'A' matrix for the additive case. 
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
T=5e3;                                      % Final time
s0 = 0.0004;                                 % strength of selection
N=1e3;                                      % Population size
Reps_s=3;                                   % No. stochastic trajectories
f_1 = [0.5; 0.25; 0.25];                    % Starting freqs stochastic replicate 1
f_2 = [0.25; 0.5; 0.25];                    % Starting freqs stochastic replicate 2
f_3 = [0.25; 0.25; 0.5];                    % Starting freqs stochastic replicate 3

% Initialize random number generator with seed for reproducibility
rng(123);

% Create matrix of fitness effects A
%A = calc_A_fluctuating(s0, n);

A = [0.010753342790922, 0.026960583349632, -0.026924388833093;
0.026960583349632, 0.006375304797180, -0.009650638297666;
-0.026924388833093, -0.009650638297666, 0.071567938794515];

% generate starting frequencies for deterministic trajectories
freqs = linspace(0,1, 51);      % generate frequencies
% Use ndgrid to create all combinations
[grid_1, grid_2, grid_3] = ndgrid(freqs, freqs, freqs);

% Reshape the grids into column vectors to form a matrix of combinations
combs = [grid_1(:), grid_2(:), grid_3(:)];

freq_combs = combs(sum(combs, 2) == 1, :);

Reps=size(freq_combs,1);        % No. of deterministic replicates

% Initialisation
Xd=zeros(n,T+1,Reps);                           % Stores deterministic trajectories
Xs=zeros(n,T+1,Reps_s);                         % Stores stochastic trajectories



% Loop for deterministic trajectories
for r=1:Reps
    Xd(1,1,r) = freq_combs(r, 1);        % Assign frequency of first allele
    Xd(2,1,r) = freq_combs(r, 2);        % Assign frequency of second allele
    Xd(3,1,r) = freq_combs(r, 3);        % Assign frequency of third allele

    % Time loop
    for k=1:T
        x=Xd(:,k,r);
        V=diag(x)-x*x';
        D=V*A*x/(1+x'*A*x);
        xp=x+D;
        xp=xp/sum(xp);                          % Forces normalisation (just to be sure!)
        Xd(:,k+1,r)=xp;
    end
end

% Loop for stochastic trajectories
for r=1:Reps_s
    % Choose starting frequencies for (first three) replicates
    if r==1
        X0= f_1; X0=X0/sum(X0);    % Forces normalisation (just to be sure!)
    elseif r==2
        X0= f_2; X0=X0/sum(X0);    % Forces normalisation (just to be sure!)
    elseif r==3
        X0= f_3; X0=X0/sum(X0);    % Forces normalisation (just to be sure!)
    else
        X0=rand(n,1); X0=X0/sum(X0);
    end

    Xs(:,1,r)=X0;

    % Time loop
    for k=1:T
        x=Xs(:,k,r);
        V=diag(x)-x*x';
        D=V*A*x/(1+x'*A*x);
        m=mnrnd(2*N,x+D)';          % Multinomial random variable (implements drift)
        xp=m/2/N;
        Xs(:,k+1,r)=xp;
    end
end



p=squeeze(abs(Xd(3,end,:)-1)<1e-2);      % Indicates trajectories where X(3) goes to 1
q=squeeze(abs(Xd(2,end,:)-0.44)<1e-2);   % Indicates trajectories where X(2) goes to 0.44

Xd1=Xd(:,:,p);                            % Set of trajectories where X(3) goes to 1               
Xd2=Xd(:,:,q);                            % Set of trajectories where X(2) goes to 0.44  

transp = 0.4;

% plotting
hold on

% deterministic loop 1
for r=1:size(Xd1,3)
    c_1 = plot3(Xd1(1,end,r),Xd1(2,end,r),Xd1(3,end,r),'r.','MarkerSize',40);     % plots final point(s) deterministic trajectories
    %hold on
    x=squeeze(Xd1(:,:,r));
    c_2 = plot3(Xd1(1,:,r),Xd1(2,:,r),Xd1(3,:,r),'linewidth',0.5);                     % plots deterministic trajectories
    set(c_2, 'Color',[0,0,1,transp]);
end

% deterministic loop 2
for r=1:size(Xd2,3)
    c_3 = plot3(Xd2(1,end,r),Xd2(2,end,r),Xd2(3,end,r),'r.','MarkerSize',40);     % plots final point(s) deterministic trajectories
    x=squeeze(Xd2(:,:,r));
    c_4 = plot3(Xd2(1,:,r),Xd2(2,:,r),Xd2(3,:,r),'linewidth',0.5);                     % plots deterministic trajectories
    set(c_4, 'Color',[1,0,0,transp]);
end

% stochastic loop
for r=1:Reps_s
    plot3(Xs(1,1,r),Xs(2,1,r),Xs(3,1,r),'b.','MarkerSize',40,'MarkerEdgeColor','k')           % plots initial point of all trajectories
    %hold on
    x=squeeze(Xs(:,:,r));
    plot3(x(1,1:T/10),x(2,1:T/10),x(3,1:T/10)+0.005,'k','linewidth',2.0)      % plots stochastic trajectories
end

box("on")
set(gca,'BoxStyle','full')
grid
axis equal
axis([0,1,0,1,0,1])

fs=35;
set(gca,'linewidth',3,'fontsize',20)
set(gca,'xtick',0:0.333:1,'xticklabel',{'0','1/3','2/3','1'})
set(gca,'ytick',0:0.333:1,'yticklabel',{'0','1/3','2/3','1'})
set(gca,'ztick',0:0.333:1,'zticklabel',{'0','1/3','2/3','1'})
set(gca,"CameraPosition",[8.4352 -2.2037 2.6731])
xlabel('$x_1$','fontsize',fs,'Interpreter','latex');
ylabel('$x_2$','fontsize',fs,'Interpreter','latex');
zlabel('$x_3$','fontsize',fs,'Interpreter','latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%
ounits = get(gcf,'Units');
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1],'Units',ounits)
orient landscape
%print -dpng Both

%saveas(gcf, 'det_and_stoch_random_regions.png');
%print("det_and_stoch_random_regions.eps",'-depsc2');
print -dpdf -bestfit -r300 det_and_stoch_random_regions
%exportgraphics(gcf,'det_and_stoch_random_regions.png','Resolution',300)
%exportgraphics(gcf,'det_and_stoch_random_regions.pdf','ContentType','image');

