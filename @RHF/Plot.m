function Plot(~, output)
ener = output.ener;
iter = output.iter;
energySet = output.energySet;

% shapes = {'s', 'o', 'v', '^', 'd', 'h', '<', '>'};
shapes = {'s', 'v', 's', 'v', 's', 'v', 's', 'v'};
colors = {'r', 'r', 'g', 'g', 'b', 'b', 'k', 'k'};

for iType = 1:length(ener)
    fprintf('%0.8f  %d \n', ener{iType}, iter{iType});
end

energyArray = [ener{:}];
iterArray = [iter{:}];
minEnergy = min(energyArray(iterArray~=100));

figure();
hold();
for iType = 1:length(ener)
    errorArray = log10(abs(energySet{iType} - minEnergy));
    errorArray(errorArray == -inf) = min(errorArray(errorArray ~= -inf));
    plot(1:length(energySet{iType}), errorArray, 'Color', colors{iType});
    scatter(1:length(energySet{iType}), errorArray, 72, shapes{iType}, colors{iType}, 'filled');
end
end
