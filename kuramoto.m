(* Kuramoto Model *)

OrderParam[graphs_, alpha_,lambda_,sigma_,K_,maxtime_,npaths_,nsteps_,seed_:1234]:= Module[{}, Export["graphs.g6",graphs];
ReadList[StringJoin[Map[ToString[#]<>" "&,{If[$OperatingSystem=="Windows","!kuramoto","!OMP_NUM_THREADS=12 ./kuramoto"],seed,npaths,nsteps,alpha,lambda,sigma,K,maxtime}]],Expression]]

glist = {RandomGraph[WattsStrogatzGraphDistribution[1000, 0.0, 5]], 
  RandomGraph[BarabasiAlbertGraphDistribution[1000, 5]], 
  RandomGraph[WattsStrogatzGraphDistribution[1000, 0.04, 5]], 
  RandomGraph[{1000, 5000}]}

maxTime = 60.0;
numPaths = 1024;
numSteps = 1024;

Print[Now];

data`Sigma1 = 
  Table[OrderParam[glist, 1.99, 0.0, sigma, 1000, maxTime, numPaths, 
    numSteps], {sigma, {0.0001, 0.01, 0.1, 0.2, 0.35, 0.5, 0.7, 1.0}}];

Print[Now];

data`Lambda1 = 
  Table[OrderParam[glist, 0.5, lambda, 0.3, 1000, maxTime, numPaths, 
    numSteps], {lambda, {0.01, 0.5, 1.5, 2.0}}];

Print[Now];

data`Lambda2 = 
  Table[OrderParam[glist, 0.5, lambda, 0.6, 1000, maxTime, numPaths, 
    numSteps], {lambda, {0.01, 0.5, 1.5, 2.0}}];

Print[Now];

data`Lambda3 = 
  Table[OrderParam[glist, 1.5, lambda, 0.3, 1000, maxTime, numPaths, 
    numSteps], {lambda, {0.01, 0.5, 1.5, 2.0}}];

Print[Now];

data`Lambda4 = 
  Table[OrderParam[glist, 1.5, lambda, 0.6, 1000, maxTime, numPaths, 
    numSteps], {lambda, {0.01, 0.5, 1.5, 2.0}}];

Print[Now];

data`Alpha1 = 
  Table[OrderParam[glist, alpha, 0.0001, 0.3, 1000, maxTime, numPaths, 
    numSteps], {alpha, {0.1, 0.5, 0.9, 0.95, 1.05, 1.1, 1.5, 1.99}}];

Print[Now];

data`Alpha2 = 
  Table[OrderParam[glist, alpha, 0.0001, 0.6, 1000, maxTime, numPaths, 
    numSteps], {alpha, {0.1, 0.5, 0.9, 0.95, 1.05, 1.1, 1.5, 1.99}}];

Print[Now];

DumpSave["results.mx", "data`"]
