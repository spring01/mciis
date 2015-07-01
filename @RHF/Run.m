function Run(obj)
[guessDensity, guessOrbital] = obj.CoreGuess();

allTypes = {'C6', 'C20', 'M6', 'M20'};
shapes = {'s', 'o', 'v', '^', 'd'};
colors = {'r', 'k', 'g', 'b', 'm'};

ener = cell(1, length(allTypes));
energySet = cell(1, length(allTypes));
iter = cell(1, length(allTypes));
for iType = 1:length(allTypes)
    [ener{iType}, energySet{iType}, iter{iType}] = obj.SCF(guessOrbital, allTypes{iType});
end

for iType = 1:length(allTypes)
    fprintf('%0.8f  %d \n', ener{iType}, iter{iType});
end

energyArray = [ener{:}];
iterArray = [iter{:}];
minEnergy = min(energyArray(iterArray~=100));

figure();
hold();
for iType = 1:length(allTypes)
    errorArray = log10(abs(energySet{iType} - minEnergy));
    errorArray(errorArray == -inf) = min(errorArray(errorArray ~= -inf));
    plot(errorArray, colors{iType});
    scatter(1:length(energySet{iType}), errorArray, 72, shapes{iType}, colors{iType}, 'filled');
end

end
