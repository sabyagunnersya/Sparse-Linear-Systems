function omegaOpt = SOROpt(rhoJacobi)

omegaOpt = 2/(1+sqrt(1-rhoJacobi^2));
fprintf('The optimum SOR factor is %f \n', omegaOpt);

end