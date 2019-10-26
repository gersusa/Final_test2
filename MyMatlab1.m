function MyMatlab1(InFile1, InFile2, InFile3, InFile4,...
    TimeLimitInSeconds, ScoringMethod, NetworkModel)
%MYMATLAB1 Version1 Script. This script generates solution1.txt
%
% Input Arguments:
% TimeLimitInSeconds and ScoringMethod are not used yet, the algorithm acts
% as it would for ScoringMethod = 0.
% InFile1: CON file
% InFile2: INL file
% InFile3: RAW file
% InFile4: ROP file


fprintf('Started\n')
startTime = tic;
%% Converting to mpc format
pool = gcp;
nworks = pool.NumWorkers;
if nworks < 2
    disp('WARNING: Parpool has less than 2 workers')
    disp('Performance will be strongly affected')
    [mpc,contingencies] = convert2mpc(InFile3,InFile4,InFile2,InFile1);
else
    [mpc,contingencies] = convert2mpc_par(InFile3,InFile4,InFile2,InFile1);
end


myfiles={'contIndex.p'};     
addAttachedFiles(gcp,myfiles); 

%% Secure-constrained Optimal Power Flow
[mpcOPF, ~, mpcOPF_or] = solveSCOPF(mpc,contingencies,...
    true,TimeLimitInSeconds, ScoringMethod, startTime);
save('mpc.mat','mpcOPF','mpcOPF_or','contingencies');
%% Create solution for short-time window
disp('Creating solution1.txt')
a=tic;
create_solution1(fixGen2Normal(gen2shunts(mpcOPF)));
toc(a)
end
