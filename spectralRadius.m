function omegaOpt = spectralRadius(GJacobi, GGS)

eJacobi = eig(GJacobi); % returns the eigenvalues 
disp('Eig Jacobi done');
eGS = eig(GGS);
disp('Eig GS done');

rhoJacobi = max(eJacobi);
rhoGS = max(eGS);

omegaOpt = SOROpt(rhoJacobi); % for calculating optimum SOR factor

fprintf('Spectral radius using Jacobi: %f \n', rhoJacobi);
fprintf('Spectral radius using GS: %f \n', rhoGS);

convRateJacobi = -log(rhoJacobi);
convRateGS = -log(rhoGS);

fprintf('Convergence rate Jacobi: %f \n', convRateJacobi);
fprintf('Convergence rate GS: %f \n', convRateGS);

disp('Spectral radius done');

end