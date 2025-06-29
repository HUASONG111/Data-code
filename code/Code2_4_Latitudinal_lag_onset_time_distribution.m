clc;
clear;close all;

%%Please change the following paths to match your own file locations.
base_path = '/Users/huasong/Desktop/data/2.Single point results/';
data_files = {'lagday_po4.csv', 'lagday_no3.csv', 'lagday_sst.csv', 'lagday_ssw.csv'};
full_paths = strcat(base_path, data_files);
titles = {'(a) PO4', '(b) NO3', '(c) SST', '(d) SSW'};
figure('Position', [100, 100, 1400, 700]);

for k = 1:length(full_paths)
    data = readmatrix(full_paths{k});
    latitudeRanges = -41:0.25:-40;

    filteredDataStruct = struct();
    
    for i = 1:length(latitudeRanges)-1
        lowerBound = latitudeRanges(i);
        upperBound = latitudeRanges(i+1);
        filteredData = data(data(:,2) >= lowerBound & data(:,2) < upperBound, :);
        fieldName = sprintf('filteredData_%.2f_%.2f', lowerBound, upperBound);
        fieldName = strrep(fieldName, '-', 'minus'); 
        fieldName = strrep(fieldName, '.', 'p');    
        filteredDataStruct.(fieldName) = filteredData;
    end
    

    fieldNames = fieldnames(filteredDataStruct);
    meanStruct = struct();
    stdStruct = struct();
    nanRatioStruct = struct();

    for i = 1:length(fieldNames)
        currentData = filteredDataStruct.(fieldNames{i});
        fourthColumn = currentData(:, 4);
        filteredValues = fourthColumn(fourthColumn >= 0); 
        meanStruct.(fieldNames{i}) = mean(filteredValues, 'omitnan');
        stdStruct.(fieldNames{i}) = std(filteredValues, 'omitnan');
        numInvalid = sum(fourthColumn == -999);
        totalValues = sum(~isnan(fourthColumn));
        nanRatioStruct.(fieldNames{i}) = numInvalid / totalValues;
    end

    xPositions = 1:length(fieldNames);
    means = cell2mat(struct2cell(meanStruct));
    stds = cell2mat(struct2cell(stdStruct));
    nanRatios = cell2mat(struct2cell(nanRatioStruct));
    
    bar_colors = {'#afd4e3','#b5d4be','#b8b9d2','#f6cbaf'};
    subplot(2, 2, k);
    yyaxis left;
    b = bar(xPositions, means, 'FaceColor', bar_colors{k}, 'EdgeColor', 'none');
    hold on;
    
    errorLineHandles = gobjects(1, length(means)); 
for i = 1:length(means)
    plot([xPositions(i), xPositions(i)], [means(i), means(i) + stds(i)], 'k-', 'LineWidth', 1);
    plot([xPositions(i) - 0.2, xPositions(i) + 0.2], [means(i) + stds(i), means(i) + stds(i)], 'k-', 'LineWidth', 1);
    errorLineHandles(i) = plot(NaN, NaN, 'k-', 'LineWidth', 1);  
end

    set(gca, 'FontSize', 10, 'FontName', 'Arial');
    xlabel('Latitude', 'FontSize', 12, 'FontName', 'Arial');
    ylabel('Mean', 'FontSize', 12, 'FontName', 'Arial');


    yyaxis right;
    l = plot(xPositions, nanRatios, 'o-', 'LineWidth', 2, 'MarkerSize', 4, ...
             'MarkerFaceColor', [227/255, 113/255, 110/255]);
    ylabel('No-lag proportion', 'FontSize', 12, 'FontName', 'Arial');


    xticks(xPositions);
    latitudeLabels = latitudeRanges(1:end-1) + diff(latitudeRanges)/2;
    xticklabels(arrayfun(@(x) sprintf('%.2f°', x), latitudeLabels, 'UniformOutput', false));

    title(titles{k}, 'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold');

    legend([b, errorLineHandles(1), l], {'Mean', 'Std', 'Ratio'}, 'Location', 'northeast');
    set(legend, 'FontSize', 12); 
    hold off;
end

drawnow;

