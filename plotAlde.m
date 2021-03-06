% SCF Done:  E(RB+HF-V5LYP) =  -153.741617377     A.U. after   16 cycles

load aldeOut.mat;

maxSCFIter = 200;
ener = output.ener;
iter = output.iter;
energySet = output.energySet;

yellow = [255 215 0] ./ 255;
darkBlue = [0 0 205] ./ 255;

% shapes = {'s', 'o', 'v', '^', 'd', 'h', '<', '>'};
shapes = {'o', 'v', 'o', 'v', 'o', 'v', 'o', 'v'};
colors = {yellow, yellow, 'k', 'k', 'g', 'g', darkBlue, darkBlue};

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
    plot(1:length(energySet{iType}), errorArray, ...
        'LineWidth', 1.5, ...
        'Color', colors{iType}, ...
        'Marker', shapes{iType}, ...
        'MarkerFaceColor', colors{iType}, ...
        'MarkerSize', 10);
end

set(gca,'FontSize', 14)

xlabel('Number of iterations', 'FontSize', 16);
ylabel('log_{10}|E_i - E_c|', 'FontSize', 16);
legend('CDIIS(6)', 'CDIIS(20)', 'MCIIS(6)', 'MCIIS(20)', 'EDIIS+CDIIS(6)', 'EDIIS+CDIIS(20)', 'EDIIS+MCIIS(6)', 'EDIIS+MCIIS(20)', ...
    'Location', 'SouthWest');

outputFileName = './graphs/alde.pdf';
print(outputFileName, '-dpdf')


