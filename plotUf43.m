

load uf43Out.mat;
% output.energySet{3} = output.energySet{3}(1:140);
% output.energySet{4} = output.energySet{4}(1:140);

maxSCFIter = 200;
ener = output.ener(3:end);
iter = output.iter(3:end);
energySet = output.energySet(3:end);

yellow = [255 215 0] ./ 255;
darkBlue = [0 0 205] ./ 255;
orange = [255 128 0] ./ 255;
darkRed = [128 0 0] ./ 255;

% shapes = {'s', 'o', 'v', '^', 'd', 'h', '<', '>'};
shapes = {'s', 'd', 's', 'd', 's', 'd'};
colors = {'g', darkBlue, orange, darkRed};

for iType = 1:length(ener)
    fprintf('%0.8f  %d \n', ener{iType}, iter{iType});
end

energyArray = [ener{:}];
iterArray = [iter{:}];
minEnergy = min(energyArray(iterArray~=maxSCFIter));

hFig = figure();
set(hFig, 'Position', [10 10 800 600])
hold();
for iType = 1:length(ener)
    errorArray = log10(abs(energySet{iType} - minEnergy));
    errorArray(errorArray == -inf) = min(errorArray(errorArray ~= -inf));
    plot(1:6:length(energySet{iType}), errorArray(1:6:end), ...
        'LineWidth', 1.5, ...
        'Color', colors{iType}, ...
        'Marker', shapes{iType}, ...
        'MarkerFaceColor', colors{iType}, ...
        'MarkerSize', 10);
end

set(gca,'FontSize', 14)

xlabel('Number of iterations', 'FontSize', 16);
ylabel('log_{10}|E_i - E_c|', 'FontSize', 16);
lgd = legend('EDIIS+CDIIS(20)', 'EDIIS+MCIIS(20)', 'ADIIS+CDIIS(20)', 'ADIIS+MCIIS(20)', ...
    'Location', 'East');

outputFileName = './graphs/uf43.pdf';
print(outputFileName, '-dpdf')


