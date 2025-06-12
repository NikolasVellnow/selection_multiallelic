% Loops over (part of) yeast genome and analyzes F for each window

clear all
close all
clc

F_results_unfiltered = load("F_results_2.mat").F_results;

% Exclude windows where optimization failed
keep_indices = ~isnan([F_results_unfiltered.cost]);
F_results = F_results_unfiltered(keep_indices);

% Exclude windows where freqs where all 0.25, F only zeros and cost = 0
F_results = F_results([F_results.cost] ~= 0);

n_wins = length(F_results);

temp_struct = struct('rep', 3, 'chr', [], ...
    'pos', [], ...
    'top_genotype', [], ...
    'type_top_genotype', [], ...
    'top_genotype_fminunc', [], ...
    'type_top_genotype_fminunc', [], ...
    'comment', []);
F_analysis_results = repmat(temp_struct, n_wins, 1);


wins_zero_cost = F_results([F_results.cost] == 0);

n_low_neg_values = 0;
n_het_top = 0;
n_het_top_fminunc = 0;
n_incongruent_top_genotype = 0;
n_incongruent_type_top_genotype = 0;

for i = 1:n_wins
    chr = F_results(i).chr;
    pos = F_results(i).pos;
    % Analysis of "final" F (either by fminunc or patternsearch)
    F = F_results(i).F;
    min_value = min(min(F));
    if min_value <= -1.0
        n_low_neg_values = n_low_neg_values +1;
    end
    
    w = convert_F_to_w(F);
    
    [max_w, linear_index] = max(w(:));
    % Convert linear index to row and column indices
    [max_row, max_col] = ind2sub(size(w), linear_index);
    F_analysis_results(i).chr = chr;
    F_analysis_results(i).pos = pos;
    genotype_label ="A_" + max_row + "A_" + max_col;
    F_analysis_results(i).top_genotype = genotype_label;
    if max_row == max_col
        type_top_genotype = "hom";
    else
        type_top_genotype = "het";
        n_het_top = n_het_top + 1;
    end
    F_analysis_results(i).type_top_genotype = type_top_genotype;

    % Analysis of F (always by fminunc)
    F_fminunc = F_results(i).F_fminunc;
    min_value = min(min(F));

    w_fminunc = convert_F_to_w(F_fminunc);
    
    [max_w_fminunc, linear_index_fminunc] = max(w_fminunc(:));
    % Convert linear index to row and column indices
    [max_row_fminunc, max_col_fminunc] = ind2sub(size(w_fminunc), linear_index_fminunc);
    genotype_label_fminunc ="A_" + max_row_fminunc + "A_" + max_col_fminunc;
    F_analysis_results(i).top_genotype_fminunc = genotype_label_fminunc;
    if max_row_fminunc == max_col_fminunc
        type_top_genotype_fminunc = "hom";
    else
        type_top_genotype_fminunc = "het";
        n_het_top_fminunc = n_het_top_fminunc + 1;
    end
    F_analysis_results(i).type_top_genotype_fminunc = type_top_genotype_fminunc;
    if genotype_label ~= genotype_label_fminunc
        n_incongruent_top_genotype = n_incongruent_top_genotype +1;
    end

    if type_top_genotype ~= type_top_genotype_fminunc
        n_incongruent_type_top_genotype = n_incongruent_type_top_genotype +1;
    end
end

percent_low_neg_values = 100*(n_low_neg_values/n_wins);
percent_het_top = 100*(n_het_top/n_wins);
fprintf("%g percent of windows have estimated fitness " + ...
    "effects smaller than -1.00\n", percent_low_neg_values);
fprintf("In %g percent of windows the most fit genotype is " + ...
    "a heterozygote\n", percent_het_top);


percent_het_top_fminunc = 100*(n_het_top_fminunc/n_wins);
fprintf("In %g percent of windows the most fit genotype is " + ...
    "a heterozygote (fminunc)\n", percent_het_top_fminunc);

% percent_incongruent_top_genotype = 100*(n_incongruent_top_genotype/n_wins);
% fprintf("In %g percent of windows the most fit genotype is " + ...
%     "incongruent between fminunc and the final F\n", percent_incongruent_top_genotype);

% percent_incongruent_type_top_genotype = 100*(n_incongruent_type_top_genotype/n_wins);
% fprintf("In %g percent of windows the type of most fit genotype is " + ...
%     "incongruent between fminunc and the final F\n", percent_incongruent_type_top_genotype);

% p1 = 100*(n_incongruent_top_genotype/n_low_neg_values);
% fprintf("In %g percent of windows with both methods the most fit genotype is " + ...
%     "incongruent between fminunc and the final F\n", p1);

%p2 = 100*(n_incongruent_type_top_genotype/n_low_neg_values);
%fprintf("In %g percent of windows with both methods the type of most fit genotype is " + ...
    %"incongruent between fminunc and the final F\n", p2);


data_table = struct2table(F_analysis_results);


% Convert the 'type_top_genotype' column to categorical
categories = categorical(data_table.chr);
top_genotype = categorical(data_table.top_genotype);

freqs = summary(top_genotype);
freqs_percent = 100*freqs.Counts/n_het_top;

% Extract the 'pos' column (numerical data)
positions = data_table.pos;

unique_top_genotypes = unique(top_genotype);

custom_colors_manual = [
    0 0 1;      % Blue
    1 0 0;      % Red
    0 0.5 1;    % Light blue
    1 0.5 0;    % Orange
    1 0.8 0;    % Yellow-orange
    0.5 0.5 1;  % Lavender/soft blue
    1 0.6 0;    % Deep Orange
    1 0.2 0.2;  % Light Red
    1 0.8 0.6;  % Peach
    0 0.3 0.8;  % Dark Blue
];

% Create custom labels with percentages
legend_labels = cell(length(unique_top_genotypes), 1);
for i = 1:length(unique_top_genotypes)
    % Find index of current genotype in freqs.Counts
    idx = find(strcmp(freqs.Categories, char(unique_top_genotypes(i))));
    
    if ~isempty(idx)
        % Create label with genotype name and percentage
        legend_labels{i} = sprintf('%s (%.1f%%)', char(unique_top_genotypes(i)), freqs_percent(idx));
    else
        % If no frequency found, just use the genotype name
        legend_labels{i} = char(unique_top_genotypes(i));
    end
end


figure;
h = gscatter(positions, categories, top_genotype, ...
    custom_colors_manual(1:length(unique_top_genotypes),:));

for i = 1:length(h)
    h(i).MarkerSize =10;
    h(i).Marker ="|";
end


xlabel('Position (bp)', 'FontSize', 20);
ylabel('Chromosome', 'FontSize', 20);
ax = gca;
ax.FontSize = 16;
legend(legend_labels{:}, 'Location', 'eastoutside');


% Adjust figure size and appearance
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7, 5]); % Set size to width=6in &; height=5in

saveas(gcf, "het_advantage_chrom_coloring.fig");
print("het_advantage_chrom_coloring.png",'-dpng','-r300');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Repeat everything with the results from fminunc
% Convert the 'type_top_genotype' column to categorical
categories = categorical(data_table.chr);
top_genotype_fminunc = categorical(data_table.top_genotype_fminunc);

freqs_fminunc = summary(top_genotype_fminunc);
freqs_percent_fminunc = 100*freqs_fminunc.Counts/n_het_top;

% Extract the 'pos' column (numerical data)
positions = data_table.pos;

unique_top_genotypes_fminunc = unique(top_genotype_fminunc);

% Create custom labels with percentages
legend_labels_fminunc = cell(length(unique_top_genotypes_fminunc), 1);
for i = 1:length(unique_top_genotypes_fminunc)
    % Find index of current genotype in freqs.Counts
    idx = find(strcmp(freqs_fminunc.Categories, char(unique_top_genotypes_fminunc(i))));
    
    if ~isempty(idx)
        % Create label with genotype name and percentage
        legend_labels_fminunc{i} = sprintf('%s (%.1f%%)', char(unique_top_genotypes_fminunc(i)), freqs_percent_fminunc(idx));
    else
        % If no frequency found, just use the genotype name
        legend_labels_fminunc{i} = char(unique_top_genotypes_fminunc(i));
    end
end


figure;
h = gscatter(positions, categories, top_genotype_fminunc, ...
    custom_colors_manual(1:length(unique_top_genotypes_fminunc),:));

for i = 1:length(h)
    h(i).MarkerSize =10;
    h(i).Marker = "|";
end

xlabel('Position (bp)', 'FontSize', 20);
ylabel('Chromosome', 'FontSize', 20);
ax = gca;
ax.FontSize = 16;
legend(h, legend_labels_fminunc{:}, 'Location', 'eastoutside');
title("fminunc");
% Adjust figure size and appearance
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7, 5]); % Set size to width=6in &; height=5in
print("het_advantage_chrom_coloring_fminunc.png",'-dpng','-r300');
