% percent change heat map
%% >> Input Parameters <<
rows = 7; %FUS
cols = 11; %RNA

%% Create Heatmap
wells = rows*cols;
percent_matrix = zeros(wells, 3);

for r = 1:rows
    for c = 1:cols
        well = (r-1)*cols + c;
        percent_matrix(well,1) = FUS(r);
        percent_matrix(well,2) = RNA(c);
        percent_matrix(well,3) = percent_increases(r,c);
    end
end
percent_tbl = array2table(percent_matrix, 'VariableNames', {'[FUS] (uM)','[RNA] (uM)','Percent Increase'});

%% Heat Map
figure;
h = heatmap(percent_tbl,'[RNA] (uM)','[FUS] (uM)','ColorVariable','Percent Increase');
h.YDisplayData = num2cell(FUS);
h.XDisplayData = num2cell(RNA);
title(strcat("Percent Increases Heatmap: ", condition))


%% Axes to Scale (Bubble Plot)
% figure;
% scatter(percent_matrix(:,2),percent_matrix(:,1),percent_matrix(:,3)*5);
% ylabel('[FUS] (uM)'), xlabel('[RNA] (uM)')
% title('Abs400 Percent Changes')