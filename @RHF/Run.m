function output = Run(obj)
[guessDensity, guessOrbital] = obj.CoreGuess();

allTypes = {'C6', 'C20', 'M6', 'M20', 'ECe6', 'ECe20', 'EMe6', 'EMe20'};
% allTypes = {'M20'};

ener = cell(1, length(allTypes));
energySet = cell(1, length(allTypes));
iter = cell(1, length(allTypes));
for iType = 1:length(allTypes)
    [ener{iType}, energySet{iType}, iter{iType}] = obj.SCF(guessOrbital, allTypes{iType});
end

output.ener = ener;
output.energySet = energySet;
output.iter = iter;

end
