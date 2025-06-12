% MultiAlleleTestTDependenceMethods.m
% The program uses 1 F0 matrix and 1 trajectory associated with this.
% It determines distance of estimated F, for different lengths of the trajectory

clear all
close all
clc

tic

rng(789);
%rng('default');

% plotting parameters
line_width = 4;
line_width_median = 8;
font_size = 50;

% Parameters
n=4;                    % No. of alleles
Tmax=50;               % Maximum length trajectory
Tset=1:1:Tmax;        % Trajectory lengths used
n_rep = 100;             % Number of replicate simulations

s_start = 0.1;         % Average fitness effect for simulated random F0s

% chose optimization method
%opt_method = "patternsearch";
opt_method = "fminunc";
%opt_method = "fmincon";

% constrained or not
%constrain = "constrained";
constrain = "unconstrained";



% Useful
nT=length(Tset);
L=n*(n+1)/2;            % No. indpendent elements in A matrix

if opt_method == "fminunc"

    % OPTIMISATION BASED ON FMINUNC
    % Set optimization options
    fminunc_options = optimoptions('fminunc', ...
            'MaxFunctionEvaluations', 100000, ...
            'TolFun', 1e-9, ...  % Function tolerance
            'TolX', 1e-9, ...    % Step tolerance
            'MaxIter', 100000, ... % Maximum number of iterations
            'Display', 'off');  % Do not display iteration information
end

if opt_method == "fmincon"

    % OPTIMISATION BASED ON FMINUNC
    % Set optimization options
    fmincon_options = optimoptions('fmincon', ...
            'MaxFunctionEvaluations', 100000, ...
            'TolFun', 1e-9, ...  % Function tolerance
            'TolX', 1e-9, ...    % Step tolerance
            'MaxIter', 100000, ... % Maximum number of iterations
            'Display', 'off');  % Do not display iteration information

warning('off', 'MATLAB:nearlySingularMatrix');
end

if opt_method == "patternsearch"
    % OPTIMISATION BASED ON PATTERNSEARCH
    ps_options = optimoptions('patternsearch', ...
            'Algorithm',"classic", ...
            'InitialMeshSize', 0.01, ...
            'MaxFunctionEvaluations', 100000, ...
            'TolFun', 1e-9, ...  % Function tolerance
            'TolX', 1e-9, ...    % Step tolerance
            'MaxIter', 100000, ... % Maximum number of iterations
            'Display', 'off');  % Do not display iteration information
end

if constrain == "constrained"
    % Contraints
    lb =-ones(L-1,1);               % Lower bound of -1  on all elements of w
    ub=Inf(L-1,1);                  % Upper bound of inf on all elements of w
end

% Initialisation
D_F=zeros(nT,n_rep);                       % D(i,r)=norm(A-F0)/norm(F0) for traj of length T=Tset(i) and replicate r
D_F_start = zeros(1,n_rep);               % To save distance before optimization start
rank_corr_F=zeros(nT,n_rep);                       
rank_corr_F_start = zeros(1,n_rep);

D_w=zeros(nT,n_rep);
D_w_start=zeros(1,n_rep);
rank_corr_w=zeros(nT,n_rep);
rank_corr_w_start=zeros(1,n_rep);


% loop over replicate simulations
for r = 1:n_rep
    fprintf('Replicate simulation: %d\n', r);
    
    X=zeros(n,Tmax+1);
    v_opt_set=zeros(L-1,nT);

    % Generate F0 matrix from a random vector
    v0=s_start*randn(L-1,1);
    F0=SMatVec([0;v0]);
    w_v0 = (v0+1)./max(v0+1);

    % How different is random F0 from matrix of zeros?
    % This is the starting point for plotting later
    %D_start(r) = norm(0.0-v0)/norm(v0);
    %D_F_start(r) = norm(0.0-v0); 
    D_F_start(r) = norm(0.0-v0/norm(v0)); 
    D_w_start(r) = norm(1.0-w_v0);
    %D_w_start(r) = norm(1.0-w_v0)/norm(w_v0);
    rank_corr_F_start(r) = 0.0;
    rank_corr_w_start(r) = 0.0;

    % Generate trajectory of length Tmax, based on F0
    aux=rand(n,1);
    X0=aux/sum(aux);                % Random initial frequency
    X(:,1)=X0;

    % loop over time points (generations)
    % Trajectory
    for t=1:Tmax
        x=X(:,t);
        V=diag(x)-x*x';
        xp=x+V*F0*x/(1+x'*F0*x);
        X(:,t+1)=xp/sum(xp);                    % Normalised at every time step (to be sure).
    end
    
    % Minimisation
    v_init=zeros(L-1,1);

    % loop over observation periods of different length
    % This simulates starting at the same point, but having
    % different number of time points in the data for estimation

    for i=1:nT
    %parfor i=1:nT
        %fprintf('Time period of length: %d\n', i);
        T=Tset(i);
        XT=X(:,1:(T+1));
        f=@(v)CostFn(v, XT);
        if opt_method == "fminunc"
            [v_opt, Cost_opt]=fminunc(f, v_init, fminunc_options);
        end
        if opt_method == "fmincon"
            [v_opt, Cost_opt]=fmincon(f, v_init, [], [], [], [], lb, ub, [], fmincon_options);
        end
        if opt_method == "patternsearch" & constrain == "unconstrained"
            [v_opt, Cost_opt]=patternsearch(f, v_init, [], [], [], [], [], [], ps_options);
        end
        if opt_method == "patternsearch" & constrain == "constrained"
            [v_opt, Cost_opt]=patternsearch(f, v_init, [], [], [], [], lb, ub, ps_options);
        end
        v_opt_set(:,i)=v_opt;
        
        %D(i,r)=norm(v_opt-v0)/norm(v0);         % A measure of distance of A and F0
        %D_F(i,r)=norm(v_opt-v0);         % A measure of distance of A and F0
        D_F(i,r)=norm(v_opt-v0)/norm(v0);
        rank_corr_F(i,r) = corr(v_opt,v0, 'Type','Spearman');
        w_v_opt = (v_opt+1)./max(v_opt+1);
        D_w(i,r) = norm(w_v_opt-w_v0);
        %D_w(i,r) = norm(w_v_opt-w_v0)/norm(w_v0);
        rank_corr_w(i,r) = corr(w_v_opt,w_v0, 'Type','Spearman');
    end
    F=SMatVec([0;v_opt]);
    disp("estimated F:")
    disp(F)

    w=SMatVec([0;w_v_opt]);
    disp("estimated w (rel fitnesses):")
    disp(w)

end

% Add distance to matrix with only zeros (null expectation)
% to the beginning of the data
D_F = [D_F_start; D_F];
rank_corr_F = [rank_corr_F_start; rank_corr_F];
D_w = [D_w_start; D_w];
rank_corr_w = [rank_corr_w_start; rank_corr_w];

D_F_median = median(D_F,2);
D_F_low_quartile = quantile(D_F,0.25,2);
D_F_high_quartile = quantile(D_F,0.75,2);

rank_corr_F_median = median(rank_corr_F,2);
rank_corr_F_low_quartile = quantile(rank_corr_F,0.25,2);
rank_corr_F_high_quartile = quantile(rank_corr_F,0.75,2);

D_w_median = median(D_w,2);
D_w_low_quartile = quantile(D_w,0.25,2);
D_w_high_quartile = quantile(D_w,0.75,2);

rank_corr_w_median = median(rank_corr_w,2);
rank_corr_w_low_quartile = quantile(rank_corr_w,0.25,2);
rank_corr_w_high_quartile = quantile(rank_corr_w,0.75,2);

Tset = [0, Tset];

toc

% Plot distance between estimated and true F matrix
figure;
plot(Tset,D_F_median,'k','LineWidth', line_width_median)
set(gca,'LineWidth',3,'fontsize', font_size)
set(gca, 'TickDir', 'out')
hold on
plot(Tset,D_F_low_quartile,':', 'Color', [0.5 0.5 0.5], 'LineWidth', line_width_median)
plot(Tset,D_F_high_quartile,':', 'Color', [0.5 0.5 0.5], 'LineWidth', line_width_median)
xlabel('length trajectory, \itt','fontsize',font_size)
ylabel('distance (F)','fontsize',font_size)

ounits = get(gcf,'Units');
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1],'Units',ounits)
orient landscape

print("FigureFDependenceT_" + opt_method + "_" + constrain ...
    + "_" + num2str(Tmax) + "_gen_distance" ...
    + ".png",'-dpng','-r300');

hold off

% Plot Spearman rank correlation between estimated and true F matrix
% entries
figure;
plot(Tset,rank_corr_F_median,'k','LineWidth', line_width_median)
set(gca,'LineWidth',3,'fontsize', font_size)
set(gca, 'TickDir', 'out')
hold on
plot(Tset,rank_corr_F_low_quartile,':', 'Color', [0.5 0.5 0.5], 'LineWidth', line_width_median)
plot(Tset,rank_corr_F_high_quartile,':', 'Color', [0.5 0.5 0.5], 'LineWidth', line_width_median)
xlabel('length trajectory, \itt','fontsize',font_size)
ylabel('Spearman rank correlation (F)','fontsize',font_size)

ounits = get(gcf,'Units');
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1],'Units',ounits)
orient landscape

print("FigureFDependenceT_" + opt_method + "_" + constrain ...
    + "_" + num2str(Tmax) + "_gen_rank_corr" ...
    + ".png",'-dpng','-r300');

hold off

% Plot distance between estimated and true rel. fitness matrix
figure;
plot(Tset,D_w_median,'k','LineWidth', line_width_median)
set(gca,'LineWidth',3,'fontsize', font_size)
set(gca, 'TickDir', 'out')
hold on
plot(Tset,D_w_low_quartile,':', 'Color', [0.5 0.5 0.5], 'LineWidth', line_width_median)
plot(Tset,D_w_high_quartile,':', 'Color', [0.5 0.5 0.5], 'LineWidth', line_width_median)
xlabel('length trajectory, \itt','fontsize',font_size)
ylabel('distance (w)','fontsize',font_size)

ounits = get(gcf,'Units');
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1],'Units',ounits)
orient landscape

print("FigureRelFitDependenceT_" + opt_method + "_" + constrain ...
    + "_" + num2str(Tmax) + "gen_distance" ...
    + ".png",'-dpng','-r300');

hold off

% Plot Spearman rank correlation between estimated and true w matrix
% entries
figure;
plot(Tset,rank_corr_w_median,'k','LineWidth', line_width_median)
set(gca,'LineWidth',3,'fontsize', font_size)
set(gca, 'TickDir', 'out')
hold on
plot(Tset,rank_corr_w_low_quartile,':', 'Color', [0.5 0.5 0.5], 'LineWidth', line_width_median)
plot(Tset,rank_corr_w_high_quartile,':', 'Color', [0.5 0.5 0.5], 'LineWidth', line_width_median)
xlabel('length trajectory, \itt','fontsize',font_size)
ylabel('Spearman rank correlation (w)','fontsize',font_size)

ounits = get(gcf,'Units');
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1],'Units',ounits)
orient landscape

print("FigureRelFitDependenceT_" + opt_method + "_" + constrain ...
    + "_" + num2str(Tmax) + "_gen_rank_corr" ...
    + ".png",'-dpng','-r300');

hold off