function fileStr = G09InputStr(info)
printOrbAlpha = [];
guessStr = [];
if(isfield(info, 'orbAlpha'))
    guessStr = ' guess=cards';
    printOrbAlpha = sprintf('%s\n', '(1e24.15)');
    for iOrb = 1:size(info.orbAlpha, 2)
        printOrbAlpha = [printOrbAlpha, sprintf('%d\n', iOrb)]; %#ok
        printOrbAlpha = [printOrbAlpha, sprintf('%24.15e\n', info.orbAlpha(:, iOrb))]; %#ok
    end
    printOrbAlpha = [printOrbAlpha, sprintf('%d\n', 0)];
end
printOrbBeta = [];
if(isfield(info, 'orbBeta'))
    for iOrb = 1:size(info.orbBeta, 2)
        printOrbBeta = [printOrbBeta, sprintf('%d\n', iOrb)]; %#ok
        printOrbBeta = [printOrbBeta, sprintf('%24.15e\n', info.orbBeta(:, iOrb))]; %#ok
    end
end

newLine = sprintf('\n');
printCommand = sprintf('%s\n', ['#p ', info.method, '/', info.basisSet , guessStr,' symmetry=none population=full scf(maxcycle=1, NoVarAcc) iop(5/33=3) iop(3/33=1) iop(5/13=1)']);
printTitle = sprintf('%s\n', 'iop(5/33=3): print fock; iop(3/33=1): print 1-e integrals; iop(5/13=1): do not terminate when scf fails');
printMol = sprintf('%d %d\n', info.chargeMult);
for iAtom = 1:size(info.cartesian, 1)
    printMol = [printMol, sprintf('%d %24.15f %24.15f %24.15f\n', info.cartesian(iAtom, :))]; %#ok
end

fileStr = [printCommand, newLine, printTitle, newLine, printMol, newLine, printOrbAlpha, printOrbBeta, newLine, newLine];
end
