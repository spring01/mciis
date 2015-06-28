function Run(obj)
[guessDensity, guessOrbital] = obj.CoreGuess();

[ener1, energySet1, iter1] = obj.SCF(guessOrbital, 'C20');
[ener2, energySet2, iter2] = obj.SCF(guessOrbital, 'M20');
[ener3, energySet3, iter3] = obj.SCF(guessOrbital, 'EC20');
[ener4, energySet4, iter4] = obj.SCF(guessOrbital, 'EM20');

fprintf('%0.8f  %d \n',ener1, iter1);
fprintf('%0.8f  %d \n',ener2, iter2);
fprintf('%0.8f  %d \n',ener3, iter3);
fprintf('%0.8f  %d \n',ener4, iter4);

energyArray = [ener1, ener2, ener3, ener4];
iterArray = [iter1, iter2, iter3, iter4];
minEnergy = min(energyArray(iterArray~=100));

figure();
hold();
errorArray1 = log10(abs(energySet1 - minEnergy));
errorArray1(errorArray1 == -inf) = min(errorArray1(errorArray1 ~= -inf));
plot(errorArray1, 'r');
scatter(1:length(energySet1), errorArray1, 72, 'square', 'r', 'filled');

errorArray2 = log10(abs(energySet2 - minEnergy));
errorArray2(errorArray2 == -inf) = min(errorArray2(errorArray2 ~= -inf));
plot(errorArray2, 'k');
scatter(1:length(energySet2), errorArray2, 72, 'o', 'k', 'filled');

errorArray3 = log10(abs(energySet3 - minEnergy));
errorArray3(errorArray3 == -inf) = min(errorArray3(errorArray3 ~= -inf));
plot(errorArray3, 'g');
scatter(1:length(energySet3), errorArray3, 72, 'v', 'g', 'filled');

errorArray4 = log10(abs(energySet4 - minEnergy));
errorArray4(errorArray4 == -inf) = min(errorArray4(errorArray4 ~= -inf));
plot(1:length(energySet4), errorArray4, 'b');
scatter(1:length(energySet4), errorArray4, 72, '^', 'b', 'filled');
end
