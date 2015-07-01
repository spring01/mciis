function fileStr = G09InputStr(chargeMult, cartesian, method, basisSet, orbital)
newLine = sprintf('\n');
printCommand = sprintf('%s\n', ['#p ', method, '/', basisSet ,' guess=cards symmetry=none population=full scf(maxcycle=1, NoVarAcc) iop(5/33=3) iop(3/33=1) iop(5/13=1)']);
printTitle = sprintf('%s\n', 'iop(5/33=3): print fock; iop(3/33=1): print 1-e integrals; iop(5/13=1): do not terminate when scf fails');
printMol = sprintf('%d %d\n', chargeMult);
for iAtom = 1:size(cartesian, 1)
    printMol = [printMol, sprintf('%d %24.15f %24.15f %24.15f\n', cartesian(iAtom, :))]; %#ok
end
printOrbital = sprintf('%s\n', '(1e24.15)');
for iOrb = 1:size(orbital, 2)
    printOrbital = [printOrbital, sprintf('%d\n', iOrb)]; %#ok
    printOrbital = [printOrbital, sprintf('%24.15e\n', orbital(:, iOrb))]; %#ok
end
fileStr = [printCommand, newLine, printTitle, newLine, printMol, newLine, printOrbital, newLine, newLine];
end
