close all; clear; clc
% written by Sophie Skanchy, version R2019b

%% Import Data
turb = readtable("WT.xlsx");

%% Prepare Data
% input parameters (<< means change the paramter)
start_well = 4; %<< 4 corresponds to column D in excel
num_rows = 7; %<< number of rows used in well plate, ex. 1 thru 6
num_cols = 11; %<< number of columns used in well plate, ex. A thru H
threshold = 10; % percent increase threshold
RNA = [8 4 2 1 .5 .25 .125 .0625 .0313 .0156 .0078]; %in uM, in order it appears on wells (descending)
FUS = [2 1.5 1 0.75 .5 .25 0]; %in uM, descending

num_wells = num_rows*num_cols;
end_col = num_wells + start_well - 1;

blank_ind = [45 46 47 56 57 58 67 68 69];
blank_col = [1 2 3];
blank_row = [5 6 7];

% process data table
turb400 = str2double(table2cell(turb(38:85, start_well:end_col)));
turb400 = turb400';
time = str2double(table2cell(turb(38:85,2)));

% create subplot titles
letters = {'RA','RB','RC','A','B','C','D','E','F','G','H'};
titles = cell(1,num_rows*num_cols);
for row = 1:num_rows
    for col = 1:num_cols
        index = (row - 1)*num_cols + col;
        titles{index} = [letters{col} num2str(row)];
    end
end



%% Plot Abs400 vs Time
max_abs = max(max(turb400));
min_abs = min(turb400(turb400>0)); %exclude blanks (0)
figure, title('Abs400 Over Time')
for i = 1:(num_cols*num_rows)
    if all(i ~= blank_ind)
        subplot(num_rows,num_cols,i)
        plot(time, turb400(i,:))
        ylim([min_abs, max_abs]) % scale to look like normalized
        title(titles{i})
    end
end

%% Percent Increases(Max/Min*100%-100%)
maxs = max(turb400');
mins = min(turb400');
percents = maxs ./ mins .* 100;
percents_matrix = zeros(num_rows, num_cols);
for row = 1:num_rows
    for col = 1:num_cols
        index = (row - 1)*num_cols + col;
        percents_matrix(row,col) = percents(index);
    end
end
percent_increases = percents_matrix - 100;

%% Percent Change 2 (tf / t0 * 100% - 100%)
% t0_abs = turb400(:,1);
% tf_abs = turb400(:,length(time));
% percents = tf_abs ./ t0_abs .* 100;
% percents_matrix = zeros(num_rows, num_cols);
% for row = 1:num_rows
%     for col = 1:num_cols
%         index = (row - 1)*num_cols + col;
%         percents_matrix(row,col) = percents(index);
%     end
% end
% percent_increases2 = percents_matrix - 100;


%% Norm Abs400
% norm_abs400 = (turb400 - min_abs)/(max_abs - min_abs);
% figure, title('Normalized Abs400 Over Time')
% for i = 1:(num_cols*num_rows)
%     subplot(num_rows,num_cols,i)
%     plot(time, norm_abs400(i,:))
%     ylim([0,1])
%     title(titles{i})
%     %display percent change
%     text(7000, .1, strcat(num2str(percent_increases(i)), "%"), "FontSize", 8)
% end

%% Percent vs FUS
figure

for i = 1:num_cols
    subplot(1,num_cols,i)
    plot(FUS, threshold*ones(1,num_rows), '-r'), hold on
    if any(i == blank_col)
        plot(FUS(1:4), percent_increases(1:4, i), '-bo')
    else
        plot(FUS, percent_increases(1:num_rows, i), '-bo')
    end
    
    ylim([0,100])
    title([num2str(RNA(i)) ' uM RNA']), xlabel('[FUS] WT (uM)')
end

subplot(1,num_cols,1), ylabel('% Increase')


%% Percent vs RNA
figure
for i = 1:num_rows
    FUS_conc = FUS(mod(i-1, num_rows) + 1);
    
    subplot(1, num_rows, i)
    plot(RNA, percent_increases(i, :), '-o'), hold on
    plot(RNA, threshold*ones(1,num_cols))
    ylim([0,100]), set(gca, 'XScale', 'log')
    title([num2str(FUS_conc) " uM FUS WT"])
    xlabel('log[RNA] (uM)')
end

subplot(1, num_rows, 1), ylabel('% Increase')

%% Run heatmap script
condition = "WT"; % for title of heatmap
run("percents_heatmap.m")
