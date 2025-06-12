% MultiAlleleTest.m

clear all
close all
clc

% Parameters
n=4;            % No. of alleles

% String of example site in genome
%ex_site_str = "rep_3_chr_C12_pos_436424";
ex_site_str = "rep_3_chr_C16_pos_457847";
%ex_site_str = "rep_3_chr_C16_pos_451847";

% read in empirical trajectory from csv file
file_str = "hap_freqs_" + ex_site_str +".csv";
trajectory_data = readtable(file_str);

% save in matrix
X = trajectory_data{:, {'DBVPG6044_freq', ...
    'DBVPG6765_freq', ...
    'Y12_freq',...
    'YPS128_freq'}}';

% generation data
g = trajectory_data{:, "time_point"}';

% Minimisation procedure
% Use fitness effects of zero as starting values
L=n*(n+1)/2;
v0=zeros(L-1,1);

% Set optimization options
fminunc_options = optimoptions('fminunc', ...
        'MaxFunctionEvaluations', 100000, ...
        'TolFun', 1e-9, ...  % Function tolerance
        'TolX', 1e-9, ...    % Step tolerance
        'MaxIter', 100000, ... % Maximum number of iterations
        'Display', 'off');  % Do not display iteration information

% Minimize CostSparseFn(v, X, g)
f=@(v)CostSparseFn(v, X, g);
[v_opt, Cost_opt]=fminunc(f, v0, fminunc_options);

F=SMatVec([0;v_opt]);
disp('Estimated F matrix (fminunc):');
disp(F);
disp('Minimum value of CostSparseFn (fminunc):');
disp(Cost_opt);

% test for non-sensical fintess values
has_value_smaller_neg_one = sum(v_opt < -1.0);

if has_value_smaller_neg_one
    ps_options = optimoptions('patternsearch', ...
                'MaxFunctionEvaluations', 100000, ...
                'TolFun', 1e-9, ...  % Function tolerance
                'TolX', 1e-9, ...    % Step tolerance
                'MaxIter', 100000, ... % Maximum number of iterations
                'Display', 'off');  % Display iteration information
    % Contraints
    lb =-ones(L-1,1);           % Lower bound of -1  on all elements of v
    ub=Inf(L-1,1);
    % Minimize CostSparseFn(v, X, g)
    f=@(v)CostSparseFn(v, X, g);
    [v_opt, Cost_opt]=patternsearch(f, v0, [], [], [], [], lb, ub, ps_options);

    F=SMatVec([0;v_opt]);
    disp('Estimated F matrix (patternsearch):');
    disp(F);
    disp('Minimum value of CostSparseFn (patternsearch):');
    disp(Cost_opt);
end


F_prime = convert_F_to_F_ij(F,3,4);
disp('Estimated F_prime with F_3,4 as reference:');
disp(F_prime);

w = convert_F_to_w(F);
disp('Estimated w matrix:');
disp(w);


%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save F matrix as a text file (might be useful later)
writematrix(F, "F_estimated_" + ex_site_str +".txt");

figure;
b = bar3(F);
title('Estimated F')


% Change observer perspective using view(azimuth, elevation)
view(-100, 20); % Example: Azimuth = -45 degrees, Elevation = 30 degrees
print("F_estimated_" + ex_site_str + ".png",'-dpng','-r300');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulating trajectories for estimated F
% Parameters
n=4;                                        % No. alleles
T=max(g);                                      % Final time
 
start_freqs = X(:, 1);               % Starting freqs from first yeast generation

% Initialisation
X_sim=zeros(n,T+1);                           % Stores deterministic trajectories
X_sim(:,1) = start_freqs;              % Assign starting freqs

% Time loop
for k=1:T
    x=X_sim(:,k);
    V=diag(x)-x*x';
    D=V*F*x/(1+x'*F*x);
    xp=x+D;
    xp=xp/sum(xp);              % Forces normalisation (just to be sure!)
    X_sim(:,k+1)=xp;
end

g_sim = 0:1:T;

font_size = 18;
font_size_ticks = 16;

% Empirical trajectories plot
figure;
hold on
plot(g, X(1,:), 'color', [0.267004, 0.004874, 0.329415], 'LineWidth', 2);   % Empirical DBVPG6044 haplotype
plot(g, X(2,:), 'color', [0.190631, 0.407061, 0.556089], 'LineWidth', 2);   % Empirical DBVPG6765 haplotype
plot(g, X(3,:), 'color', [0.20803, 0.718701, 0.472873], 'LineWidth', 2);   % Empirical Y12 haplotype
plot(g, X(4,:), 'color', [0.993248, 0.906157, 0.143936], 'LineWidth', 2);   % Empirical YPS128 haplotype

title('Yeast experimental evolution', 'FontSize', font_size);
xlabel('Generations', 'FontSize', font_size);
ylabel('Haplotype frequency', 'FontSize', font_size);
ax = gca;
ax.FontSize = font_size_ticks;
ylim([0.0,1.0]);
xlim([0, 18]);

legend({'DBVPG6044','DBVPG6765','Y12','YPS128'}, 'Location', 'northwest');

% Adjust figure size and appearance
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 6, 5]); % Set size to width=6in &; height=5in
print("trajectories_empirical_" + ex_site_str + ".png",'-dpng','-r300');


% Modelled trajectories plot
figure;
hold on
plot(g_sim, X_sim(1,:), 'color', [0.267004, 0.004874, 0.329415], 'LineWidth', 2);   % Simulated DBVPG6044 haplotype
plot(g_sim, X_sim(2,:), 'color', [0.190631, 0.407061, 0.556089], 'LineWidth', 2);   % Simulated DBVPG6765 haplotype
plot(g_sim, X_sim(3,:), 'color', [0.20803, 0.718701, 0.472873], 'LineWidth', 2);   % Simulated Y12 haplotype
plot(g_sim, X_sim(4,:), 'color', [0.993248, 0.906157, 0.143936], 'LineWidth', 2);   % Simulated YPS128 haplotype

title('Model trajectories', 'FontSize', font_size);
xlabel('Generations', 'FontSize', font_size);
ylabel('Haplotype frequency', 'FontSize', font_size);
ax = gca;
ax.FontSize = font_size_ticks;
ylim([0.0,1.0]);
xlim([0, 18]);

legend({'DBVPG6044','DBVPG6765','Y12','YPS128'}, 'Location', 'northwest');

% Adjust figure size and appearance
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 6, 5]); % Set size to width=6in &; height=5in
print("trajectories_model_" + ex_site_str + ".png",'-dpng','-r300');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% predicting trajectories for estimated A for 22 more generations
% Parameters
n=4;                                        % No. alleles
T=22;                                      % Final time
 
start_freqs = X(:, 17);               % Starting freqs from yeast generation 540

% Initialisation
X_sim=zeros(n,T+1);                           % Stores deterministic trajectories
X_sim(:,1) = start_freqs;              % Assign starting freqs

% Time loop
for k=1:T
    x=X_sim(:,k);
    V=diag(x)-x*x';
    D=V*F*x/(1+x'*F*x);
    xp=x+D;
    xp=xp/sum(xp);              % Forces normalisation (just to be sure!)
    X_sim(:,k+1)=xp;
end

g_sim = max(g):1:(max(g)+T);

figure;
hold on

% first plot empirical trajectories of experimental evolution
plot(g, X(1,:), ...
    'color', [0.267004, 0.004874, 0.329415], ...
    'LineWidth', 2);   % DBVPG6044 haplotype
plot(g, X(2,:), ...
    'color', [0.190631, 0.407061, 0.556089], ...
    'LineWidth', 2);   % DBVPG6765 haplotype
plot(g, X(3,:), ...
    'color', [0.20803, 0.718701, 0.472873], ...
    'LineWidth', 2);   % Y12 haplotype
plot(g, X(4,:), ...
    'color', [0.993248, 0.906157, 0.143936], ...
    'LineWidth', 2);   % YPS128 haplotype

% plot predicted trajectories for follwoing 22 generations
plot(g_sim, X_sim(1,:), ...
    'color', [0.267004, 0.004874, 0.329415], ...
    'LineWidth', 2, ...
    'LineStyle',':');   % Simulated DBVPG6044 haplotype
plot(g_sim, X_sim(2,:), ...
    'color', [0.190631, 0.407061, 0.556089], ...
    'LineWidth', 2, ...
    'LineStyle',':');   % Simulated DBVPG6765 haplotype
plot(g_sim, X_sim(3,:), ...
    'color', [0.20803, 0.718701, 0.472873], ...
    'LineWidth', 2, ...
    'LineStyle',':');   % Simulated Y12 haplotype
plot(g_sim, X_sim(4,:), ...
    'color', [0.993248, 0.906157, 0.143936], ...
    'LineWidth', 2, ...
    'LineStyle',':');   % Simulated YPS128 haplotype

title('Prediction (for generation > 18)', 'FontSize', font_size);
xlabel('Generations', 'FontSize', font_size);
ylabel('Haplotype frequency', 'FontSize', font_size);
ax = gca;
ax.FontSize = font_size_ticks;
ylim([0.0,1.0]);
xlim([0, max(g)+T]);

legend({'DBVPG6044','DBVPG6765','Y12','YPS128'}, 'Location', 'northwest');
grid on;


% Adjust figure size and appearance
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 6, 5]); % Set size to width=6in &; height=5in
print("trajectories_predicted_" + ex_site_str + ".png",'-dpng','-r300');



% analyse estimated matrix
disp('"Total heterosis" for estimated matrix w:')
w_homs = diag(w);
mean_w_homs = mean(w_homs);
fprintf('Average fitness of homozygotes: %.3f\n', mean_w_homs);

mask = triu(true(size(w)), 1);
w_hets = w(mask);
mean_w_hets = mean(w_hets);
fprintf('Average fitness of heterozygotes: %.3f\n', mean_w_hets);

fprintf('Heterozygote fitness is %.3f times as high as homozygote fitness\n', mean_w_hets/mean_w_homs);


% analyse estimated matrix
disp('"Pairwise heterosis" for estimated matrix w:')

ph_vector = [];

for allel_1 = 1:n
    for allel_2 = 1:n
        if allel_1 > allel_2
            w_het = w(allel_1,allel_2);
            w_homs = (w(allel_1,allel_1)+w(allel_2,allel_2))/2;
            ph = w_het/w_homs;
            fprintf('heterosis for genotype A%dA%d: %.3f\n',allel_1, allel_2, ph);
            ph_vector = [ph_vector, ph];
        end
    end
end
