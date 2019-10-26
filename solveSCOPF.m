function [ mpcOPF, pfs, mpcOPF_or, result ] = solveSCOPF (mpc, contingencies,...
    verbose, TimeLimitInSeconds, ScoringMethod, startTime )
%SOLVESCOPF Solve security-constrained OPF
%   [MPCOPF, PFS, MPCOPF_OR, RESULT] = SOLVESCOPF(MPC, CONTINGENCIES,
%   VERBOSE)
%
%   MPCOPF is the OPF result of the last iteration before exiting the loop.
%   If successful, this is the solution to the SCOPF problem.
%
%   PFS is an array of MatPower cases that contains each post-contingency
%   MatPower system. Is empty if OPF does not converge.
%
%   MPCOPF_OR is the OPF result of the first iteration. In other words,
%   this is the solution of the non-security-constrained OPF problem.
%
%   RESULT is the RESULT output argument of the RUNALLCONS function in the
%   last iteration before exiting the loop. If OPF does not converge, it
%   only contains two vectors with -1 values (voltViol and noConvs).
%
%   MPC is a MatPower case. It must be the result of CONVERT2MPC, after
%   extending the OPF, i.e. it must be the output of the
%   TOGGLE_SOFTLIMS_EXT(MPC,'on') function.
%
%   CONTINGENCIES must be the CONTINGENCIES struct returned by the
%   CONVERT2MPC function.
%
%   VERBOSE is an optional input. If true, information is printed each
%   iteration. False by default.
%
%   About current version: 10

%      The new pfs struct is in external notation, with the fields:
%      pfs(i).bus     [Index_bus VM1 VA1 VMAX VMIN BS1]
%      pfs(i).branch  [From_bus To_bus Pf Qf Pt Qt]
%      pfs(i).gen     [Gen_bus PG QG VG]

%   VoltPenalty is now  a output of runALLCons, it was necessary to cut the
%   pfs fields respect the 9 version 


%   Current version: 8
%
%   About this version:
%   Fixed persistent bug when contingencies.gen or contingencies.branch are
%   empty.
%   Modified verbose behavior.
%
%   About version 7:
%   Fixed bug when contingencies.gen or contingenies.branch is empty
%
%   About version 6:
%   Modified stage 'overloads' of finite state machine. Now generator
%   limits are also modified in generator contingencies.
%   Fixed bug that appears when there is no active swing bus in the system
%   at the start.
%
%   About version 5:
%   Included the call to EXTEND_OPF within SOLVESCOPF. Adapted to code to
%   new OPF extensions for area spinning reserves and reactive shunts.
%
%   About version 4:
%   Changed solver to IPOPT.
%   Now using switched shunts converted to generators.
%   Modified voltages and noConvs states in FSM according to new
%   formulation (May 29 2019)
%   Extended verbose output.
%   Coded possibility to use one iteration OPF as seed of the next, but
%   left it commented out for now.
%
%   About version 3:
%   Added finite state machine to advance through algorithm stages
%   (overloads and voltages).
%
%   About version 2.1:
%   Fixed bug. Output arguments PFS and RESULT correspond to best solution,
%   i.e. to the MPCOPF returned.
%
%   About version 2:
%   mpcOPF returned is now the best solution attained throughout all 
%   iterations. This is done because there is no guarantee that the last
%   iteration is the best solution.
%
%   Copyright (c) 2019, Gers USA
%   by 
%   Tomas Valencia tvalencia@gersusa.com
%   Daniel Agudelo-Martinez daniel.agudelo@gers.com.co
%   Dario Arango dario.arango@gers.com.co
%   Camilo Acosta camilo.acosta@gers.com.co
% 
%% Validate inputs
if nargin<3
    % By default, do not print iteration information
    verbose = false;
end

% Start time
if nargin<6
    startTime = tic;
end
if ~isa(startTime,'uint64') || ~isscalar(startTime)
    error('Invalid start time')
end

if nargin<5
    % Div 1 by default
    ScoringMethod = 1;
end
if nargin<4
    % Div 1 by default
    TimeLimitInSeconds = 60000e2;
end


%% Constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
%Constants PFS now
VM1=2;
VA1=3;
VMAX1=4;
VMIN1=5;


%% Parameters
% Tolerance for all hard limits
TOL = 1e-4;
% Checking contingencies structure
if ~isempty(contingencies.branch)
   conK = keys(contingencies.branch);
else
    conK = [];
end
if ~isempty(contingencies.gen)
    conG = keys(contingencies.gen);
else
    conG = [];
end
%%-- Parameters dependent on size of system
% Seed: automatic for Div2, mpc case for Div1
% Contingency screening trigger
selCont = false;
rtCont=1;
if ScoringMethod == 2 || ScoringMethod == 4
    %Div 2 or Div 4: -network size- criteria
    if size(mpc.bus,1) >= 15e3
        selCont = true;
    end
else
    %Div 1 or Div 2: -number of contingencies- criteria
    if length(conK)+length(conG) > 1000
        selCont = true;
    end    
end
% Time margin & estimation for first iteration
% Initial state of FSM
% Improvement threshold in overloads stage of FSM
if size(mpc.bus,1) < 2e3
    tMargin = 7;
    initState = 'overloads';
    improvThreshold = -1;
    selCont = false;
elseif size(mpc.bus,1) >= 2e3 && size(mpc.bus,1) < 3e3
    tMargin = 10;
    initState = 'voltages';
    improvThreshold = -20;
    selCont = false;
elseif size(mpc.bus,1) >= 3e3 && size(mpc.bus,1) < 5e3
    tMargin = 10;
    initState = 'voltages';
    improvThreshold = -100;
    if selCont
        % Contingency selection equation
        rtCont=(-0.05)*((length(conK)+length(conG))/...
            size(mpc.bus,1))+0.725;
    end
elseif size(mpc.bus,1)>=5e3 && size(mpc.bus,1)<8e3
    tMargin = 15;
    initState = 'voltages';
    improvThreshold = -200;
    if selCont
        % Contingency selection equation
        rtCont=(-0.07)*((length(conK)+length(conG))/...
            size(mpc.bus,1))+0.725;
    end
else
    tMargin = 20;
    initState = 'voltages';
    improvThreshold = -500;
    if selCont
        if size(mpc.bus,1) <= 15e3 ||...
                (ScoringMethod == 2 || ScoringMethod == 4)
            % Contingency selection equation
            rtCont=(-0.1)*((length(conK)+length(conG))/...
                size(mpc.bus,1))+0.725;
        else % buses>15k & Div1
            % Contingency selection equation
            if size(mpc.bus,1) < 30e3
                rtCont=0.1;
                %secPerCont = 0.2;
            else
                rtCont=0.05;
                %secPerCont = 0.3;
            end
        end
    end
end
rtCL=round(rtCont,2);% Rate of branch CONs evaluated
rtCG=round(rtCont,2);% Rate of gen CONs evaluated


%% Initialization
% MPC_CP will contain the OPF problem that will be modified each iteration
% but is not overwritten with the OPF solution
mpc_cp = mpc;

%-- Miscellaneous variables
tOPF1 = TimeLimitInSeconds - toc(startTime) - 5*tMargin;
seedFdbk = true; %Use current case as seed for next OPF
result.voltViol = [-1 -1];
result.noConvs = [-1 -1];
ctgPen = [];
voltPenalty = []; 
pfs = [];
it = 0; % Iteration counter
stage = initState;
getout = false; % Break loop
tTotSCOPF = tic;
improv = -Inf(1,2);
lastScore = Inf;
vvPerf = Inf(1,2);
t_post_RAC=30;
tOPFfactor = 1.3;
tSelCont0=30;
limitTimeRAC = TimeLimitInSeconds- tMargin;
%--- Parameters for parallel scripts
pool=gcp;
pool.IdleTimeout = 120;
disp(strcat('ParPool IdleTime: ',num2str(pool.IdleTimeout),'minutes'))
numWorks=pool.NumWorkers;
if numWorks < 7
    disp('WARNING: Parpool has less than 7 workers')
    disp('Performance will be strongly affected')
end
comb_num = min(numWorks, 30); % 45 available combinations for ipopt.options
% Loading different set of ipopt options
[OPT, accept_index]= loadSeeds(comb_num, tOPF1);

if verbose
    fprintf(['Iteration,Cost,Penalty,TotCost,VoltPenalty,VoltViol',...
        ',NoConvs,Improv1,Improv2,',...
        'tOPF1,s1,tOPF2,s2,tOPF3,s3,tOPF4,s4,tRAC,tPen,tSave,Stage,',...
        'tFSM1,tFSM2,tTot,rlvntCont\n'])    
end

%---- Convert switched shunts to generators and extend OPF
mpc_cp=swshunt2gen(mpc_cp);
mpc_cp = extend_opf(mpc_cp,'on',contingencies);


%% Main Loop
while it < 40
    
    tTot = tic;
    %-- Set variables only for first iteration if not enough time to finish
    if it==0
        mpcOPF = mpc_cp;
        mpcOPF.f = 1e15-1; % Dummy cost
        mpcOPF_best = mpcOPF;
        mpcOPF_or = mpcOPF;
        pfs_best = pfs;
        result_best = result;
        ctgPen_best = ctgPen;
        voltPenalty_best = voltPenalty;
    end
    % Check remaining time
    if toc(startTime) > TimeLimitInSeconds-(tOPF1+tMargin)
        fprintf('Timeout\n')
        break
    end
    % Check there is a swing bus
    if ~any(mpc_cp.bus(:,BUS_TYPE)==REF)
        ag = find(mpc_cp.gen(:,GEN_STATUS)~=0);
        [~,refg] = max(mpc_cp.gen(ag,PMAX));
        mpc_cp.bus(mpc_cp.bus(:,BUS_I)==mpc_cp.gen(ag(refg(1)),GEN_BUS),...
            BUS_TYPE)=REF;
    end
    
    %%---------------------------------------------------------------------
    % 1. Run parOPF
    %%---------------------------------------------------------------------
    % OPF parallel version
    % Sending tasks to workers
    tOPF1 = tic;
    F(1:comb_num)=parallel.FevalFuture;
    for k=1:comb_num
        F(k)=parfeval(@myrunOPF,1,mpc_cp,OPT(k));
    end
    %-- Checking worker status
    runningW=[];
    successFlag = false;
    % Either one OPF was successful or all OPFs were unsuccessful
    while ~successFlag
        runningW=strcmp({F(1:comb_num).State},'running');
        if ~all(runningW)
            try
                Fout=[F(~runningW).OutputArguments];
                Fout_mat=cell2mat(Fout(:));
                success = [Fout_mat(:).success];
                if any(success ~= 0) || sum(~success)==comb_num
                    successFlag = true;
                end
            catch ME
                disp('Error in parOPF execution:')
                disp(ME.message)
                [mpcOPF_best,mpcOPF_or] = update_finish(mpc_cp,mpcOPF,...
                    mpcOPF_best,mpcOPF_or,it);
                getout = true;
                delete(gcp('nocreate'))
                break
            end
        end
    end
    if getout
        break
    end
    %-- Waiting for better solutions if required
    if any(accept_index(~runningW))
        disp('Waiting for better solutions')
        pause(5);
    else
        pause(0.5);
    end
    %-- Cancel remaining tasks on other workers
    cancel(F(1:comb_num));
    
    %-- Retrieve objective function from parOPF
    try
        cancelledW_idx=cellfun(@(x) ~isempty(x), {F(1:comb_num).Error});
        OPFfuture=fetchOutputs(F(~cancelledW_idx));
        cost=[OPFfuture(:).f];
        [~,ind_opf]=min(cost,[],2);
        mpcOPF=OPFfuture(ind_opf);
    catch ME
        disp(strcat('Error in parOPF execution: '),ME.cause)
        delete(gcp('nocreate'))
        break
    end
    s1 = mpcOPF.success;
    tOPF1 = toc(tOPF1)
    %-- Set variables only for first iteration
    if it==0
        mpcOPF_or = mpcOPF;
        mpcOPF_best = mpcOPF;
        itEnd = tOPF1;
    else
        t_post_RAC=itEnd-tRAC;
    end
    %-- Check remaining time
    if toc(startTime)>TimeLimitInSeconds-(itEnd+tMargin)
        fprintf('Timeout after OPF1\n')
        [mpcOPF_best,mpcOPF_or] = update_finish(mpc_cp,mpcOPF,mpcOPF_best...
            ,mpcOPF_or,it);
        break
    end
    delete(F(:))
    %clearvars F
    
    
    itStart=tic;
    % Legacy vars
    tOPF2 = -1;
    s2 = -1;
    tOPF3 = -1;
    s3 = -1;
    tOPF4 = -1;
    s4 = -1;
    
    %-- Deactivate OPF extensions for contingency evaluation. Otherwise a
    % lot of errors occur
    if toggle_softlims_ext(mpc_cp,'status')
        mpcOPF = toggle_softlims_ext(mpcOPF,'off');
    end
    if toggle_qflow(mpc_cp,'status')
        mpcOPF = toggle_qflow(mpcOPF,'off');
    end
    if toggle_qshuntsgen(mpc_cp,'status')
        mpcOPF = toggle_qshuntsgen(mpcOPF,'off');
    end
    if toggle_areaspinreserve(mpc_cp,'status')
        mpcOPF = toggle_areaspinreserve(mpcOPF,'off');
    end
    
    if it==0
        % Keep original OPF for limit check
        mpcOPF_or = mpcOPF;
    end
    
    % If OPF does not converge, exit loop
    if ~mpcOPF.success
        if verbose
            fprintf('OPF NoConv, out\n')
        end
        %if ~exist('mpcOPF_best','var')
        if (mpcOPF.f < mpcOPF_best.f) && it<2
            mpcOPF_best = mpcOPF;
        end
        pfs_best = pfs;
        result_best = result;
        mpc_cp_best = mpc_cp;
        voltPenalty_best = voltPenalty;
        %end
        break
    end
    
    
    
    
    %%---------------------------------------------------------------------
    % 2. Contingency evaluation
    %%---------------------------------------------------------------------
    
    % If not enough time left, quit
    if toc(startTime)>TimeLimitInSeconds-(itEnd+tMargin)
        if verbose
            fprintf('Timeout before selCont\n')
        end
        if (mpcOPF.f <= 10*mpcOPF_best.f) && it<2
            mpcOPF_best = mpcOPF;
        end
        break
    end
    
    %%--Contingency Selection
    tSelCont = 0;
    
    if selCont
        
        if it == 0
            fprintf('selCont activated\n')
            fprintf('Branch rate: %5.2f\n',rtCL);
            fprintf('Gens rate: %5.2f\n',rtCG);
            tSelCont=tSelCont0;
            contingencies_m=contingencies;    
        end
        %-- Rate estimation
        Nb=length(mpcOPF.bus(:,1));
        tpc=time_per_cont(Nb);
        N_tot=length(contingencies.branch)+length(contingencies.gen);
        %Define contingences rate
        limitTimeRAC=TimeLimitInSeconds-toc(startTime)-tOPFfactor*tOPF1...
            -t_post_RAC-tSelCont;
        rtCL = min(rtCL,(limitTimeRAC/(tpc*N_tot))*(numWorks/72))
        rtCG=rtCL;
        tSelCont = tic;
        % Serial or parallel execution
        if numWorks < 7 || size(mpcOPF.bus,1)>15e3
            disp('selectContingencies NO_PAR')
            conIndex = selectContingencies_nf(rtCL,rtCG,mpcOPF,...
                contingencies_m,'p');
        else
            disp('selectContingencies PAR')
            conIndex = selectContingencies_nf_par3(rtCL,rtCG,mpcOPF,...
                contingencies_m,'p');
        end
        contingencies = editContingencies(conIndex,contingencies_m);
        % Redefinition of conK, conG
        if isfield(contingencies,'branch') && ~isempty(contingencies.branch)
            conK = keys(contingencies.branch);
        else
            contingencies.branch=[];
            conK = [];
        end
        if isfield(contingencies,'gen') && ~isempty(contingencies.gen)
            conG = keys(contingencies.gen);
        else
            contingencies.gen=[];
            conG = [];
        end
       
        tSelCont = toc(tSelCont)
    end
    
    % If not enough time left, quit
    toc(startTime)
    if toc(startTime)>TimeLimitInSeconds-(itEnd-tSelCont+tMargin)
        if verbose
            fprintf('Timeout before runALLCONS\n')
        end
        if (mpcOPF.f <= 10*mpcOPF_best.f) && it < 2
            mpcOPF_best = mpcOPF;
        end
        break
    end
    
    % Run all contingencies on current base case and check
       
    if isempty(contingencies.branch) && isempty(contingencies.gen)
    fprintf('Timeout after selCont\n')
        [mpcOPF_best,mpcOPF_or] = update_finish(mpc_cp,mpcOPF,mpcOPF_best,mpcOPF_or,it);
        break
    end 
    
    tRAC = tic;
    
     %limitTimeRAC = TimeLimitInSeconds-(itEnd-tSelCont+tMargin);
   [result, contingencies ,pfs, voltPenalty]= runAllCONS(mpcOPF,contingencies,...
        mpcOPF_or,'AC',0,limitTimeRAC);
    tRAC = toc(tRAC)
    length(contingencies.branch)+length(contingencies.gen)
    tRAC/(length(contingencies.branch)+length(contingencies.gen))
    % % %         if selCont
    % % %             selCont_timeFactor = rtCL / tRAC;
    % % %             remTime = TimeLimitInSeconds - toc(startTime);
    % % %             rtCL = round(min(rtCL, selCont_timeFactor *...
    % % %                 (remTime-tFSM1-tFSM2-tOPF1)), 2);
    % % %             rtCG = rtCL;
    % % %         end
    % Contingency penalty: overload + power mismatch
    tPen = tic;
    ctgPen=0.5.*sum(evalPenalty(cat(1,result.sigma.sigma_c.sigmaTh{:})))./...
        (length(contingencies.branch)+length(contingencies.gen))+...
        0.5.*sum(result.sigma.sigmaMis(:))./...
        (length(contingencies.branch)+length(contingencies.gen));
    tPen = toc(tPen);
    
    
    
    
    %%---------------------------------------------------------------------
    % 3. Save best solution
    %%---------------------------------------------------------------------
    tSave = tic;
    if it==0
        % Initialize best solution
        mpcOPF_best = mpcOPF;
        ctgPen_best = ctgPen;
        result_best = result;
        pfs_best = pfs;
        mpc_cp_best = mpc_cp;
        voltPenalty_best = voltPenalty;
        contingencies_best=contingencies;
    end
    
    % Update best solution (without taking voltages into account)
    if (mpcOPF.f + ctgPen < mpcOPF_best.f + ctgPen_best) ...
            && sum(1*result.voltViol)<=sum(1*result_best.voltViol)...
            && sum(1*result.noConvs)<=sum(1*result_best.noConvs)
        mpcOPF_best = mpcOPF;
        ctgPen_best = ctgPen;
        result_best = result;
        pfs_best = pfs;
        mpc_cp_best = mpc_cp;
        voltPenalty_best = voltPenalty;
        contingencies_best=contingencies;
    end
    tSave = toc(tSave);
    
    % Evaluate improvement of current solution
    improv = circshift(improv,1,2);
    improv(1)=mpcOPF.f+ctgPen-lastScore;
    lastScore = mpcOPF.f + ctgPen;
    if isequal(initState,'voltages') && mpcOPF.f > (1e3*mpcOPF_best.f)
        if verbose
            fprintf('Over restricted Network, out\n')
        end
        break
    end
    
    if verbose
        fprintf(['%d,%6.2f,%7.2f,%7.2f,%7.2f,%d,%d,%7.2f',...
            ',%7.2f,',...
            '%4.0f,%d,%4.0f,%d,%4.0f,%d,%4.0f,%d,',...
            '%4.0f,%4.0f,%4.0f,%s,'],...
            it,mpcOPF.f,ctgPen,mpcOPF.f+ctgPen,voltPenalty,...
            sum(result.voltViol),sum(result.noConvs),...
            improv(1),improv(2),...
            tOPF1,s1,tOPF2,s2,tOPF3,s3,tOPF4,s4,...
            tRAC,tPen,tSave,stage)
    end
    
    
    
    
    %%---------------------------------------------------------------------
    % 4. Finite state machine
    %%---------------------------------------------------------------------
    % If not enough time left, quit
    if toc(startTime)>TimeLimitInSeconds-(itEnd-tSelCont...
            -tRAC-tPen-tSave+tMargin)
        if verbose
            fprintf('Timeout before FSM\n')
        end
        break
    end
    tFSM1 = tic;
    %%-- Transitions
    switch stage
        case 'overloads'
            % If no more overloads or no more progress advance to voltage
            % stage
            if isempty(result.sigma.sigmaTh.summary) ||...
                    all(improv>improvThreshold)
                
                mpc_cp = mpc_cp_best;
                mpcOPF = mpcOPF_best;
                result = result_best;
                ctgPen = ctgPen_best;
                voltPenalty = voltPenalty_best;
                contingencies=contingencies_best;
                clearvars pfs
                pfs=pfs_best;
                if any(result.voltViol)
                    stage = 'voltages';
                elseif any(result.noConvs)
                    stage = 'noConvs';
                else
                    getout = true; % Break loop
                    stage = 'exit';
                end
            else
                % To do: address OL and voltViol simultaneously
                if any(result.voltViol) && ~isequal(initState,'overloads')
                    stage = 'voltages';
                elseif any(result.noConvs)
                    stage = 'noConvs';
                end
            end
        case 'voltages'
            % Update best solution taking voltages into account
            if mpcOPF.f + ctgPen + voltPenalty <...
                    mpcOPF_best.f + ctgPen_best + voltPenalty_best
                mpcOPF_best = mpcOPF;
                ctgPen_best = ctgPen;
                result_best = result;
                voltPenalty_best = voltPenalty;
                clearvars pfs_best
                pfs_best = pfs;
                mpc_cp_best = mpc_cp;
                contingencies_best=contingencies;
            end
            
            % If no more voltage violations, end
            if ~any(result.voltViol)
                mpcOPF_best = mpcOPF;
                ctgPen_best = ctgPen;
                result_best = result;
                clearvars pfs_best
                pfs_best = pfs;
                mpc_cp_best = mpc_cp;
                voltPenalty_best = voltPenalty;
                contingencies_best=contingencies;
                if any(result.noConvs)
                    stage = 'noConvs';
                else
                    if ctgPen > 0.05*mpcOPF.f
                        % NOTE: This 5% threshold could be reviewed
                        stage = 'overloads';
                        improv = -Inf(1,2);
                    else
                        stage = 'exit';
                        getout = true; % Break loop
                    end
                end
            else
                % To do: address OL and voltViol simultaneously
                %                 if any(result.noConvs)
                %                     stage = 'noConvs';
                %                 elseif any(result.noConvs)
                %                     stage = 'noConvs';
                %                 end
            end
        case 'noConvs'
            % If no more NoConvs, end
            if ~any(result.noConvs)
                mpcOPF_best = mpcOPF;
                ctgPen_best = ctgPen;
                result_best = result;
                clearvars pfs_best
                pfs_best = pfs;
                mpc_cp_best = mpc_cp;
                voltPenalty_best = voltPenalty;
                contingencies_best=contingencies;
                
                if any(result.voltViol)
                    stage = 'voltages';
                else
                    if ctgPen > 0.05*mpcOPF.f
                        % NOTE: This 5% threshold could be reviewed
                        stage = 'overloads';
                        improv = -Inf(1,2);
                    else
                        stage = 'exit';
                        getout = true; % Break loop
                    end
                end
                
            end
    end
    tFSM1 = toc(tFSM1);
    tFSM2 = tic;
    %%- Modify OPF of next iteration, depending on current state
    
    switch stage
        case 'overloads'
            %%-- Update line ratings to compute next OPF
            summary = result.sigma.sigmaTh.summary;
            getFlow = @(x) max(sqrt(x(1).^2+x(2).^2),sqrt(x(3).^2+x(4).^2));
            for ol = 1:size(result.sigma.sigmaTh.summary,1)
                c = summary(ol,3);
                %-- Overloaded branch
                s_ol = getFlow(mpcOPF.branch(summary(ol,1),[PF,QF,PT,QT]));
                
                if c<=length(contingencies.branch)
                    %-- Branch contingency
                    br = mpcOPF.indexMap.branch.psse2mpc(contingencies.branch(conK{c}));
                    s_br = getFlow(mpcOPF.branch(br,[PF,QF,PT,QT]));
                    delta_s = result.sigma.sigmaTh.val(c);
                    %%--
                    if s_ol>s_br
                        br_redef = summary(ol,1);
                        s_redef = s_ol;
                    else
                        br_redef = br;
                        s_redef = s_br;
                    end
                    mpc_cp.branch(br_redef,RATE_A) = min(...
                        s_redef,mpc_cp.branch(br_redef,RATE_A)) -...
                        delta_s*mpcOPF.baseMVA;
                else
                    %-- Generator contingency
                    gen = mpcOPF.indexMap.gen.psse2mpc(contingencies.gen(conG{c-...
                        length(contingencies.branch)}));
                    pg = mpcOPF.gen(gen,PG);
                    qg = mpcOPF.gen(gen,QG);
                    sg = sqrt(pg^2+qg^2);
                    
                    delta_s = result.sigma.sigmaTh.val(c)*mpcOPF.baseMVA;
                    
                    delta_pg = delta_s*pg/sg;
                    delta_qg = delta_s*qg/sg;
                    
                    %-- Limit capacity of largest element between generator
                    %and OL branch
                    if s_ol>sg ||...
                            ((mpcOPF.branch(summary(ol,1),RATE_A)-s_ol)<TOL...
                            && summary(ol,2)<0.5)
                        br_redef = summary(ol,1);
                        s_redef = s_ol;
                        mpc_cp.branch(br_redef,RATE_A) = min(...
                            mpc_cp.branch(br_redef,RATE_A),s_redef) -...
                            delta_s;
                    else
                        mpc_cp.gen(gen,PMAX) = max(...
                            mpc_cp.gen(gen,PMIN),...
                            pg - delta_pg);
                        if delta_qg>0
                            mpc_cp.gen(gen,QMAX) = max(qg - delta_qg,...
                                mpc_cp.gen(gen,QMIN));
                        else
                            mpc_cp.gen(gen,QMIN) = min(qg - delta_qg,...
                                mpc_cp.gen(gen,QMAX));
                        end
                    end
                end
            end
            
        case 'voltages'
            if ~isempty(contingencies.branch)
                conBV = values(contingencies.branch);
                % Vector for storing branch CONs that cause
                % non-marginal voltage violations
                lv_b = result.voltViol(1:length(conBV));
                % Vector for storing gen CONs that cause
                % non-marginal voltage violations
                lv_g = result.voltViol(length(conBV)+1:end);
                
                TOL_2 = 25E-3; % Tolerance for defining marginal voltage violations
                outMat_max = []; % Matrix for storing new bus VMAX limits
                outMat_min = []; % Matrix for storing new bus VMIN limits
                tao=1;
                
                % Function that computes power factor of (P,Q) pair
                powfact = @(x) cos(atan2(x(:,2),x(:,1)));
                
                for i=find(result.voltViol)'
                    %%-- Look for marginal voltage violations
                    %mpcInt = int2ext(pfs(i));
                    mpcInt =pfs{i};
                    busV_max = mpcInt.bus(:,VM1)-mpc.bus(:,VMAX)>TOL & ...
                        mpcInt.bus(:,VM1)-mpc.bus(:,VMAX)<TOL_2 & ...
                        abs(mpcOPF.bus(:,VM)-mpc.bus(:,VMAX))<TOL_2;
                    busV_max_dV = mpcOPF.bus(busV_max,VM) - ...
                        (mpcInt.bus(busV_max,VM1)-mpc.bus(busV_max,VMAX));
                    outMat_max = [outMat_max; [find(busV_max), busV_max_dV]];
                    busV_min = mpc.bus(:,VMIN)-mpcInt.bus(:,VM1)>TOL & ...
                        mpc.bus(:,VMIN)-mpcInt.bus(:,VM1)<TOL_2 & ...
                        abs(mpcOPF.bus(:,VM)-mpc.bus(:,VMIN))<TOL_2;
                    busV_min_dV = mpcOPF.bus(busV_min,VM) + ...
                        (mpc.bus(busV_min,VMIN)-mpcInt.bus(busV_min,VM1));
                    outMat_min = [outMat_min; [find(busV_min), busV_min_dV]];
                    
                    % If all violations in this CON were small, do not include
                    % in il nor ig
                    if all(busV_max==...
                            (mpcInt.bus(:,VM1)-mpc.bus(:,VMAX)>TOL)) && ...
                            all(busV_min==...
                            (mpc.bus(:,VMIN)-mpcInt.bus(:,VM1)>TOL))
                        if i>length(conBV)
                            lv_g(i-length(conBV))=false;
                        elseif i<=length(conBV)
                            lv_b(i)=false;
                        end
                    end
                    
                    %%- Look for branches with high power factor
                    if i>length(conBV)
                        % If gen contingency, skip
                        continue
                    end
                    % If line has high PowFact, modify voltage limit at base
                    % case for line ends that present voltage violations
                    il = mpc.indexMap.branch.psse2mpc(conBV{i});
                    pfF = powfact(mpcOPF.branch(il,[PF,QF]));
                    pfT = powfact(mpcOPF.branch(il,[PT,QT]));
                    f_bus = mpcInt.bus(:,BUS_I)==mpcOPF.branch(il,F_BUS);
                    ovF = mpcInt.bus(f_bus,VM1) - mpc.bus(f_bus,VMAX);
                    uvF = mpc.bus(f_bus,VMIN) - mpcInt.bus(f_bus,VM1);
                    t_bus = mpcInt.bus(:,BUS_I)==mpcOPF.branch(il,T_BUS);
                    ovT = mpcInt.bus(t_bus,VM1) - mpc.bus(t_bus,VMAX);
                    uvT = mpc.bus(t_bus,VMIN) - mpcInt.bus(t_bus,VM1);
                    %                     vvBus = mpcInt.bus(mpcInt.bus(:,VM)>mpc.bus(:,VMAX) | ...
                    %                         mpcInt.bus(:,VM)<mpc.bus(:,VMIN),BUS_I);
                    vvBus = mpcInt.bus(mpcInt.bus(:,VM1)>mpcInt.bus(:,VMAX1)+TOL | ...
                        mpcInt.bus(:,VM1)<mpcInt.bus(:,VMIN1)-TOL,BUS_I);
                    if abs(pfF)>0.95 && ismember(mpc.bus(f_bus,BUS_I),vvBus)
                        if ovF>TOL
                            outMat_max = [outMat_max;...
                                [find(f_bus),1+0.8*(mpcOPF.bus(f_bus,VM)-1)]];
                        elseif uvF>TOL
                            outMat_min = [outMat_min;...
                                [find(f_bus),1+0.8*(mpcOPF.bus(f_bus,VM)-1)]];
                        end
                        %                         lv_b(i)=false;
                        tao=1.05;
                    end
                    if abs(pfT)>0.95 && ismember(mpc.bus(t_bus, BUS_I),vvBus)
                        if ovT>TOL
                            outMat_max = [outMat_max;...
                                [find(t_bus),1+0.8*(mpcOPF.bus(t_bus,VM)-1)]];
                        elseif uvT>TOL
                            outMat_min = [outMat_min;...
                                [find(t_bus),1+0.8*(mpcOPF.bus(t_bus,VM)-1)]];
                        end
                        %                         lv_b(i)=false;
                        tao=1.05;
                    end
                    
                    %                     if (abs(pfF)||abs(pfT)>0.95) && sum(mpcOPF.branch(il,[QF,QT]))<0
                    %                         mpc_cp.bus(ismember(mpc_cp.bus(:,BUS_I),vvBus),VMAX) =...
                    %                             1+0.8*(mpcOPF.bus(ismember(mpcOPF.bus(:,BUS_I),vvBus),VM)-1);
                    %                         lv_b(i)=false;
                    %                     end
                    %
                end
                % Transfer info from outMat_max and outMat_min to mpc_cp
                for busV = unique(outMat_max(:,1))'
                    mpc_cp.bus(busV,VMAX)=...
                        min(outMat_max(outMat_max(:,1)==busV,2));
                end
                for busV = unique(outMat_min(:,1))'
                    mpc_cp.bus(busV,VMIN)=...
                        max(outMat_min(outMat_min(:,1)==busV,2));
                end
                
                
                %%-- Modify OPF for next iteration for contingencies causing
                %%non-marginal and non-over-.95PF voltage issues
                
                % Get outage branches in non-marginal voltage violation CONs
                il = cellfun(@(x) mpc.indexMap.branch.psse2mpc(x),conBV(...
                    logical(lv_b)))';
                
                % Add reactive limit to selected outage branches
                if ~isempty(il)
                    vvPerf = circshift(vvPerf,1,2);
                    vvPerf(1) = sum(il)+sum(vvBus);
                    if isequal(vvPerf(1),vvPerf(2)) && ~any(isinf(vvPerf))
                        QBrLim = QBrLim*0.9;
                    else
                        QBrLim=0.9;
                    end
                    % Set limit at 90% of old limit
                    mpc_cp.branch(il,RATE_B) = min((QBrLim*tao).*min(...
                        abs(mpcOPF.branch(il,[QF,QT])),[],2),...
                        mpc_cp.branch(il,RATE_B));
                    if exist('qflow','var') && isfield(qflow,'il')
                        qflow.il = unique([il;qflow.il]);
                    else
                        qflow.il = il;
                    end
                    mpc_cp.qflow = qflow;
                    % Extend OPF to add reactive flow constraints
                    if ~toggle_qflow(mpc_cp,'status')
                        mpc_cp = toggle_qflow(mpc_cp,'on');
                    end
                end
            end
            
            % Get outage generators in non-marginal voltage violation CONs
            if ~isempty(contingencies.gen)
                conGV = values(contingencies.gen);
                ig = cellfun(@(x) mpc.indexMap.gen.psse2mpc(x),conGV(...
                    logical(lv_g)))';
                % Set reactive limit of outage generators at narrower values
                if ~isempty(ig)
                    % Set limit at 80% of old QG
                    vvBus=mpcInt.bus(mpcInt.bus(:,VM1) > ...
                        mpcInt.bus(:,VMAX1)+TOL,BUS_I);
                    if ~isempty(vvBus)
                        if any(ismember(vvBus, mpcOPF.bus(ismember(...
                                mpcOPF.bus(:,BUS_I),mpcOPF.gen(ig,...
                                GEN_BUS)),BUS_I)))
                            posQ = mpcOPF.gen(ig,QG)>0;
                            %                     vvBus_g=mpcInt.bus(mpcInt.bus(:,VM)>mpc.bus(:,VMAX) | ...
                            %                         mpcInt.bus(:,VM)<mpc.bus(:,VMIN),BUS_I);
                            %                     if ismember(mpc.bus(ismember(mpc.bus(:,BUS_I),mpc.gen(ig(posQ),GEN_BUS)),BUS_AREA),...
                            %                         mpcInt.bus(vvBus_g,BUS_AREA))
                            mpc_cp.gen(ig(posQ),QMAX) = ...
                                0.8*mpcOPF.gen(ig(posQ),QG);
                            mpc_cp.gen(ig(~posQ),QMIN) = ...
                                0.8*mpcOPF.gen(ig(~posQ),QG);
                        else
                            mpc_cp.bus(ismember(mpc_cp.bus(:,BUS_I),...
                                vvBus),VMAX) = 1+0.8*(mpcOPF.bus(...
                                ismember(mpcOPF.bus(:,BUS_I),vvBus),VM)-1);
                            %                     else
                            %                         posP = mpcOPF.gen(ig,PG)>0;
                            %                         mpc_cp.gen(ig(posP),PMAX) = 0.8*mpcOPF.gen(ig(posP),PG);
                            %                         mpc_cp.gen(ig(~posP),PMIN) = 0.8*mpcOPF.gen(ig(~posP),PG);
                            %                     end
                        end
                    else
                        %%--undervoltage
                    end
                end
            end
            
        case 'noConvs'
            %%-- Modify OPF of next iteration to prevent noConvs
            
            % For branch CON, reduce flow through branch in base case to
            % 90% of old value
            if ~isempty(contingencies.branch)
                conBV = values(contingencies.branch);
                il = cellfun(@(x) mpc.indexMap.branch.psse2mpc(x),conBV(...
                    logical(result.noConvs(1:length(conBV)))))';
                getFlow = @(x) max([sqrt(x(:,1).^2+x(:,2).^2),...
                    sqrt(x(:,3).^2+x(:,4).^2)],[],2);
                s_il = getFlow(mpcOPF.branch(il,[PF,QF,PT,QT]));
                mpc_cp.branch(il,RATE_A) = 0.9.*s_il;
            end
            
            % For gen CONs, reduce PMAX to 75% of old PG value
            % NOTE: Reduction of QG currently commented out, might be
            % enabled again if deemed necessary
            if ~isempty(contingencies.gen)
                conGV = values(contingencies.gen);
                ig = cellfun(@(x) mpc.indexMap.gen.psse2mpc(x),conGV(...
                    logical(result.noConvs(length(conBV)+1:end))))';
                if ~isempty(ig)
                    posP = mpcOPF.gen(ig,PG)>0;
                    mpc_cp.gen(ig(posP),PMAX) = ...
                        max(0.75*mpcOPF.gen(ig(posP),PG),...
                        mpc_cp.gen(ig(posP),PMIN));
                    mpc_cp.gen(ig(~posP),PMIN) = ...
                        max(0.75*mpcOPF.gen(ig(~posP),PG),...
                        mpc_cp.gen(ig(~posP),PMIN));
                end
            end
    end
    tFSM2 = toc(tFSM2);
    tTot = toc(tTot);
    if verbose
        fprintf('%4.0f,%4.0f,%4.0f\n',tFSM1,tFSM2,tTot)
    end
    
    
    %-- Transfer info from mpc_cp to mpcOPF for using result of current
    %iteration as seed of next iteration.
    %-- Update seed for next OPF
    if seedFdbk
        % Retrieve data from current OPF
        VAL=mpcOPF.var.val;
        x0=[VAL.Va;
            VAL.Vm;
            VAL.Pg;
            VAL.Qg;
            VAL.y;
            VAL.s_rate_a_lin;
            VAL.s_rate_a_xfmr;
            VAL.s_p_unbal_plus;
            VAL.s_p_unbal_minus;
            VAL.s_q_unbal_plus;
            VAL.s_q_unbal_minus;
            VAL.s_ar
            ];
        % Assign current data to mpc_cp
        mpc_cp.x0=x0;
        mpc_cp.lamb0_e=[ mpcOPF.nle.lambda.Pmis ; mpcOPF.nle.lambda.Qmis ] ;
        mpc_cp.lamb0_i=[ mpcOPF.nli.mu.softSf_lin ; mpcOPF.nli.mu.softSt_lin;...
            mpcOPF.nli.mu.softSf_xfmr; mpcOPF.nli.mu.softSt_xfmr;...
            mpcOPF.nli.mu.Qsh] ;
        % Adapt current seed if qflow_toogle was triggered
        if isfield(mpcOPF.nli.mu,'Qf')
            mpc_cp.lamb0_i= [mpc_cp.lamb0_i;mpcOPF.nli.mu.Qf];
            mpc_cp.lamb0_i= [mpc_cp.lamb0_i;mpcOPF.nli.mu.Qt];
        end
    end
    
    clearvars mpcOPF OPF
    
    %-- Update iteration counter
    it = it+1;
    itEnd = toc(itStart);
    
    %-- Break loop if applicable
    if getout
        break
    end
end
%% Function Output
% Return best solution attained w/o swshunts converted back
% mpcOPF is not converted back so that Code2 has access to all info
mpcOPF = mpcOPF_best;%fixGen2Normal(gen2shunts(mpcOPF_best));

% Convert back to external indexing for creating output files and convert
% switched shunts back
% NOTE: This must be commented out in submitted version of the code
% % % if ~isempty(pfs)
% % %     for i=1:length(pfs)
% % %         pfs(i) = fixGen2Normal(gen2shunts(int2ext(pfs(i))));
% % %     end
% % % end

if verbose
    tTotSCOPF = toc(tTotSCOPF);
    pfs = pfs_best;
    result = result_best;
    ctgPen = ctgPen_best;
    voltPenalty = voltPenalty_best;

        fprintf('%d,%6.2f,%7.2f,%7.2f,%7.2f,%d,%d,%4.0f\n',...
            -1,mpcOPF.f,ctgPen,mpcOPF.f+ctgPen,voltPenalty,...
            sum(result.voltViol),sum(result.noConvs),...
            tTotSCOPF)
end
end

%% ========================================================================
%%-- Nested functions

%-- myrunOPF function
function mpcOPF = myrunOPF( mpc,OPT)
%-- Define constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA,...
    BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
TOL=1e-4;
%--Input parameters
if  OPT.Vm_ones_flag
    mpc.bus(:,VM)=1;
end
seed=OPT.seed;
optionsipopt=OPT.optionsipopt;
%-- OPF function
mpcOPF = runopf(mpc, mpoption('verbose',0,'out.all',0,...
    'opf.violation',TOL,'opf.ac.solver','ipopt',...
    'opf.start',seed,'model','AC', 'ipopt.opts', optionsipopt));
end



%-- Timeout Function
function [mpcOPF_best,mpcOPF_or] = update_finish(mpc_cp, mpcOPF,...
    mpcOPF_best, mpcOPF_or, it)

% Update best solution (without taking voltages into account)
if isfield(mpcOPF,'f') && (mpcOPF.f <= 10*mpcOPF_best.f) && (it < 2)
    mpcOPF_best = mpcOPF;
end
%-- Deactivate OPF extensions for contingency evaluation. Otherwise a
% lot of errors occur
if toggle_softlims_ext(mpc_cp,'status')
    mpcOPF_best = toggle_softlims_ext(mpcOPF_best,'off');
    mpcOPF_or = toggle_softlims_ext(mpcOPF_or,'off');
end
if toggle_qflow(mpc_cp,'status')
    mpcOPF_best = toggle_qflow(mpcOPF_best,'off');
    mpcOPF_or = toggle_qflow(mpcOPF_or,'off');
end
if toggle_qshuntsgen(mpc_cp,'status')
    mpcOPF_best = toggle_qshuntsgen(mpcOPF_best,'off');
    mpcOPF_or = toggle_qshuntsgen(mpcOPF_or,'off');
end
if toggle_areaspinreserve(mpc_cp,'status')
    mpcOPF_best = toggle_areaspinreserve(mpcOPF_best,'off');
    mpcOPF_or = toggle_areaspinreserve(mpcOPF_or,'off');
end
end
function alfa=time_per_cont(Nb)

%retorna el tiempo por contingencia que se tarda un for  en funcion del 
%numero de nodos.
Nb=Nb+1;
NBx=  [4918   8733  10000  11615 19402 30000];
alfas=[0.023  0.028  0.037  0.041  0.178 0.363];

if (Nb<NBx(1)) %Nodos menores
%Se supone la misma pendiente que entre Nb(1) y Nb(2)
m=(alfas(2)-alfas(1))/(NBx(2)-NBx(1));
alfa=alfas(1)-m*(NBx(1)-Nb);
elseif (Nb>NBx(end)) %Nodos mayores
%Se supone la misma pendiente que entre Nb(end) y Nb(end-1)
m=(alfas(end)-alfas(end-1))/(NBx(end)-NBx(end-1));
alfa=alfas(end)+m*(Nb-NBx(end));
else
[~,min]=find(Nb>NBx,1, 'last');
[~,max]=find(Nb<=NBx,1);
m=(alfas(max)-alfas(min))/(NBx(max)-NBx(min));
alfa=alfas(min)+m*(Nb-NBx(min));
end
alfa
end