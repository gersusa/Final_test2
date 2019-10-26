function [conIndex] = selectContingencies_nf_par2(rtCL,rtCG,mpcOPF,contingencies,option)
%SELECTCONTINGENCIES select a subset of contingencies
%   [CONINDEX] = SELECTCONTINGENCIES(RTCL, RTCG,
%   MPCOPF,CONTINGENCIES,OPTION)
%
%   CONINDEX contains the indices of selected contingencies
%
%   RTCL is the selection rate in p.u. for contingencies.branch
%
%   RTCG is the selection rate in p.u. for contingencies.gen
%
%   MPCOPF is the OPF result.
%
%   CONTINGENCIES must be the CONTINGENCIES struct returned by the
%   CONVERT2MPC function.
%
%   OPTION allows to select different strategies for contingencies
%   selection
%
%   Current version: 6
%
%   About this version:
%   Parallelization is used on sorted lists computation.
%
%   About version 5:
%   'For' loops were changed by nested functions. contIndex function now
%   converts from element index to contingency index.
%
%   About version 4:
%   Sorted lists of elements are created using additional cirteria:
%   Pbr(t,f), Qbr(t,f), Pgn, Qgn
%
%   About version 3:
%   The empty option (' ') is included in main switch. This option allows
%   to select 100% of contigencies.
%
%   About version 2:
%   Check for empty fields in contingencies strcuture and retrieve
%   contingencies length
%
%   Copyright (c) 2019, Gers USA
%   by Daniel Agudelo-Martinez daniel.agudelo@gers.com.co

%% Constants
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
BASE_KV=10;
BUS_I=1;
%%

% Check input variables
if ~isempty(contingencies.branch)
    conKeysB = keys(contingencies.branch);
    lconKB = length(conKeysB);
else
    lconKB = 0 ;
end
if ~isempty(contingencies.gen)
    conKeysG = keys(contingencies.gen);
    lconKG = length(conKeysG);
else
    lconKG = 0;
end

%Different criteria for contingencies selection
%conIndex_l=0;
%Selecting ordered contigencies
%Branch contingencies

    conBV=values(contingencies.branch);
    getFlow = @(x) max([sqrt(x(:,1).^2+x(:,2).^2),...
        sqrt(x(:,3).^2+x(:,4).^2)],[],2);
    powfact = @(x) cos(atan2(x(:,2),x(:,1)));
    sList_BR_out=[];
    mpcOPF_branch=mpcOPF.branch;
    mpcOPF_indexMap_branch=mpcOPF.indexMap.branch;
    mpcOPF_bus = mpcOPF.bus;
    pool = gcp;
    nworks = pool.NumWorkers
    numlabs

    
    
    parfor_time = tic;
    parfor i=1:7
        %             sList_BRs=[];
        %             sList_BRq=[];
        %             sList_BRqD=[];
        %             sList_BRv=[];
        %             sList_BRpf=[];
        %             sList_BRpfD=[];
        %             sList_BRconn=[];
        % sBR_ix_S = [];
        % sBR_ix_Q
        % sBR_ix_QD
        
        
        
        if i==1
            tic
            % Sorted lines according to real power flow (from)
            %             getFlow = @(x) max([sqrt(x(:,1).^2+x(:,2).^2),...
            %                 sqrt(x(:,3).^2+x(:,4).^2)],[],2);
            
            %             BrS=[find(mpcOPF.branch(:,1)) getFlow(mpcOPF.branch(:,[PF,QF,PT,QT]))];
            %             XfS=sortrows([BrS(mpcOPF.isxfmr,1) BrS(mpcOPF.isxfmr,2)],-2);
            %             sList_xs= contIndex( mpcOPF.indexMap.branch, XfS(:,1), conBV );
            %             LnS=sortrows([BrS(~mpcOPF.isxfmr,1) BrS(~mpcOPF.isxfmr,2)],-2);
            %             sList_ls= contIndex( mpcOPF.indexMap.branch, LnS(:,1), conBV );
            [~, sBR_ix_S]=...
                sort(getFlow(mpcOPF_branch(:,[PF,QF,PT,QT])),'descend');
            %sList_BRs= contIndex( mpcOPF_indexMap_branch, sBR_ix_S, conBV );
            sList_BR_out(:,i) = contIndex( mpcOPF_indexMap_branch, sBR_ix_S, conBV );
            toc
        elseif i==5
            %             ------------------------------------------------------------------
            
            %             % Sorted lines according to real power flow (to)
            %             [~, sLines_ix_pt]=sortrows(abs(mpcOPF.branch(mpcOPF.branch(:,11)>0,:)),-16);
            %             conIndex_l_pt= contIndex( mpcOPF.indexMap.branch, sLines_ix_pt, conBV );
            % 			conIndex_l_pt=conIndex_l_pt(1:round(length(conIndex_l_pt)*rtCL));
            
            tic
            % Sorted lines according to reactive power (from)
            %[~, sLines_ix_qf]=sort(abs(mpcOPF.branch(:,QF)),'descend');
            %sList_lqf= contIndex( mpcOPF.indexMap.branch, sLines_ix_qf, conBV );
            %sList_l=sList_l(1:round(length(sList_lq)*rtCL));
            %                                     -------------------------------------------------
            %                                     BrQ=[find(mpcOPF.branch(:,1)) max(mpcOPF.branch(:,[QF,QT]),[],2)];
            %                                     XfQ=sortrows([BrQ(mpcOPF.isxfmr,1) BrQ(mpcOPF.isxfmr,2)],-2);
            %                                     sList_xq= contIndex( mpcOPF.indexMap.branch, XfQ(:,1), conBV );
            %                                     LnQ=sortrows([BrQ(~mpcOPF.isxfmr,1) BrQ(~mpcOPF.isxfmr,2)],-2);
            %                                     sList_lq= contIndex( mpcOPF.indexMap.branch, LnQ(:,1), conBV );
            [~, sBR_ix_Q]=...
                sort(max(abs(mpcOPF_branch(:,[QF,QT])),[],2),'descend');
            %sList_BRq= contIndex( mpcOPF_indexMap_branch, sBR_ix_Q, conBV );
            sList_BR_out(:,i) = contIndex( mpcOPF_indexMap_branch, sBR_ix_Q, conBV );
            toc
        elseif i==3
            tic
            % Sorted lines according to Q diff
            %[~, sLines_ix_qt]=sort(abs(mpcOPF.branch(:,QT)),'descend');
            %sList_lqt= contIndex( mpcOPF.indexMap.branch, sLines_ix_qt, conBV );
            %sList_l=sList_l(1:round(length(sList_l)*rtCL));
            [~, sBR_ix_QD]=...
                sort(abs(diff(abs(mpcOPF_branch(:,[QF,QT])),1,2)),'descend');
            %sList_BRqD= contIndex( mpcOPF_indexMap_branch, sBR_ix_QD, conBV );
            sList_BR_out(:,i) = contIndex( mpcOPF_indexMap_branch, sBR_ix_QD, conBV );
            toc
        elseif i==4
            tic
            %                         % Sorted lines according to power rate
            %                         [~, sLines_ix_r]=sort(abs(mpcOPF.branch(:,RATE_A)),'descend');
            %                         sList_l= contIndex( mpcOPF.indexMap.branch, sLines_ix_r, conBV );
            %                         sList_l=sList_l(1:round(length(sList_l)*rtCL));
            % Sorted lines according to voltage
            %                                     sLines_ix_v=[];
            %                                     sBus_ix_v=unique([mpcOPF.bus(mpcOPF.bus(:,10)>100,1); mpcOPF.bus(mpcOPF.bus(:,10)>100,2)]);%from&to
            %                                     for ii=1:length(sBus_ix_v)
            %                                     sLines_ix_v=unique([sLines_ix_v;find(mpcOPF.branch(mpcOPF.branch(:,11)==1,1)==sBus_ix_v(ii))]);
            %                                     end
            %                                     conIndex_l_v= contIndex( mpcOPF.indexMap.branch, sLines_ix_v, conBV );
            
            % % %                                 for i=1:size(mpcOPF.branch,1)
            % % %                                     voltFBus(i,2)=mpcOPF.bus(find(mpcOPF.branch(i,2)==mpcOPF.bus(:,1)),10);
            % % %                                 end
            %---Tests-----------------------
            %                                 voltFBus=arrayfun(@(x) mpcOPF.bus(find(x==mpcOPF.bus(:,BUS_I),1),BASE_KV), mpcOPF.branch(:,F_BUS));
            %                                 voltTBus=arrayfun(@(x) mpcOPF.bus(find(x==mpcOPF.bus(:,BUS_I),1),BASE_KV), mpcOPF.branch(:,T_BUS));
            [voltFBus, voltTBus] = voltBus_function(mpcOPF_branch, mpcOPF_bus, F_BUS, T_BUS, BUS_I, BASE_KV);
            %                                 -------Working only for netowrks having
            %                                 less or equal branches than buses-------
            %                                 [~,sBR_ix_V]=sort(max([mpcOPF.bus(mpcOPF.branch(:, F_BUS),BASE_KV),...
            %                                    mpcOPF.bus(mpcOPF.branch(:,T_BUS),BASE_KV)],[],2),'descend');
            %                                 ----------------------------------------
            [~,sBR_ix_V]=sort(max([voltFBus,voltTBus],[],2),'descend');
            %                                 [~,sBR_ix_V]=sort(100*(max([mpcOPF.bus(mpcOPF.branch(:,F_BUS),BASE_KV), ...
            %                                     mpcOPF.bus(mpcOPF.branch(:,T_BUS),BASE_KV)],[],2)<300),'ascend');
            %sList_BRv= contIndex( mpcOPF_indexMap_branch, sBR_ix_V, conBV );
            sList_BR_out(:,i) = contIndex( mpcOPF_indexMap_branch, sBR_ix_V, conBV );
            toc
        elseif i==2
            tic
            %                                  Sorted lines according to power factor
            pfT=abs(powfact(mpcOPF_branch(:,[PT,QT])));
            pfF=abs(powfact(mpcOPF_branch(:,[PF,QF])));
            [~,sBR_ix_pf]=sort(min([pfF,pfT],[],2),'ascend');%pf value
            %sList_BRpf= contIndex( mpcOPF_indexMap_branch, sBR_ix_pf, conBV );
            sList_BR_out(:,i) = contIndex( mpcOPF_indexMap_branch, sBR_ix_pf, conBV );
            toc
        elseif i==6
            tic
            pfT=abs(powfact(mpcOPF_branch(:,[PT,QT])));
            pfF=abs(powfact(mpcOPF_branch(:,[PF,QF])));
            [~,sBR_ix_pfD]=sort(abs(pfF-pfT),'descend');%pf differences
            %sList_BRpfD= contIndex( mpcOPF_indexMap_branch, sBR_ix_pfD, conBV );
            sList_BR_out(:,i) = contIndex( mpcOPF_indexMap_branch, sBR_ix_pfD, conBV );
            toc
        elseif i==7
            tic
            %                                 Sorted lines according to connectivity
            % % %                                 for i=1:size(mpcOPF.branch,1)
            % % %                                     connBRF(i,1)=sum(mpcOPF.branch(:,1)==mpcOPF.branch(i,1));
            % % %                                 end
            %------Tests--------
            %                                 connBRF=arrayfun(@(x) sum(mpcOPF.branch(:,F_BUS)==x), mpcOPF.branch(:,F_BUS));
            %                                 connBRT=arrayfun(@(x) sum(mpcOPF.branch(:,T_BUS)==x), mpcOPF.branch(:,T_BUS));
            [connBRF, connBRT] = connBR_function(mpcOPF_branch, F_BUS, T_BUS);
            [~,sBR_ix_conn]=sort(max([connBRF connBRT],[],2),'ascend');%connectivity
            
            %sList_BRconn= contIndex( mpcOPF_indexMap_branch, sBR_ix_conn, conBV );
            sList_BR_out(:,i) = contIndex( mpcOPF_indexMap_branch, sBR_ix_conn, conBV );
            toc
            %                                 otherwise
            %                         %sList_l=[];
        end
    end
    
    
    %             [~,rList_BRs]=sort([sList_BRs{:}]);
    %             [~,rList_BRq]=sort([sList_BRq{:}]);
    %             [~,rList_BRqD]=sort([sList_BRqD{:}]);
    %             [~,rList_BRv]=sort([sList_BRv{:}]);
    %             [~,rList_BRpf]=sort([sList_BRpf{:}]);
    %             [~,rList_BRpfD]=sort([sList_BRpfD{:}]);
    %             [~,rList_BRconn]=sort([sList_BRconn{:}]);
    parfor_time= toc(parfor_time)
    [~,rList_BRs]=sort(sList_BR_out(:,1));
    [~,rList_BRq]=sort(sList_BR_out(:,5));
    [~,rList_BRqD]=sort(sList_BR_out(:,3));
    [~,rList_BRv]=sort(sList_BR_out(:,4));
    [~,rList_BRpf]=sort(sList_BR_out(:,2));
    [~,rList_BRpfD]=sort(sList_BR_out(:,6));
    [~,rList_BRconn]=sort(sList_BR_out(:,7));
    
    
    %             aux=[rList_BRs, rList_BRq, rList_BRpf, rList_BRpfD, rList_BRconn, rList_BRv];
    w=[0.01,0.01,1,2,1, 1, 2]';
    wMean=w/sum(w);
    % % % %             [~,aux_list2]=sort(mean([...
    % % % %                 rList_BRs...
    % % % %                 ,rList_BRq...
    % % % %                 ,rList_BRpf...
    % % % %                 ,rList_BRpfD...
    % % % %                 ,rList_BRconn...
    % % % %                 rList_BRv...
    % % % %                 ],2),'ascend');
    [~,aux_list2]=sort([...
        rList_BRs...
        ,rList_BRq...
        ,rList_BRpf...
        ,rList_BRpfD...
        ,rList_BRconn...
        ,rList_BRv...
        ,rList_BRqD...
        ]*wMean,'ascend');
    %List_xfmr=find(mpcOPF.isxfmr);
    %sList_l=List_xfmr;
    %sList_l=unique([sList_l;aux_list2]);
    sList_l=aux_list2;
    
    conIndex_l=sList_l(1:round(length(sList_l)*rtCL));



%Gen contingencies
if lconKG>0
    conGV=values(contingencies.gen);
    
    getS = @(x) sqrt(x(:,1).^2+x(:,2).^2);
    sGen = getS(mpcOPF.gen(:,[QG,PG]));
    [~,sGEN_ix_s] = sort(sGen,'descend');
    i_r = ismember(sGEN_ix_s,mpcOPF.reactgen);
    sList_GENs= contIndex(mpcOPF.indexMap.gen, sGEN_ix_s(~i_r), conGV )+lconKB;
    
    %             Sorted gens according to power factor
    powfact = @(x) cos(atan2(x(:,2),x(:,1)));
    pfGen=abs(powfact(mpcOPF.gen(:,[PG,QG])));
    [~,sGEN_ix_pf]=sort(pfGen,'ascend');
    i_r = ismember(sGEN_ix_pf,mpcOPF.reactgen);
    sList_GENpf= contIndex(mpcOPF.indexMap.gen, sGEN_ix_pf(~i_r), conGV )+lconKB;
    
    
    
    [~,rList_GENs]=sort(sList_GENs);
    [~,rList_GENpf]=sort(sList_GENpf);
    
    [~,aux_list3]=sort(min([...
        rList_GENs...
        ,rList_GENpf...
        ],[],2),'ascend');
    
    %sList_g=sList_gs;
    sList_g=aux_list3+lconKB;
    conIndex_g=sList_g(1:round(length(sList_g)*rtCG));
    
else
    conIndex_g=0;
end



%Branch and Gen contingency indices
conIndex=[conIndex_l;...
    conIndex_g];
end
%% ========================================================================
%-- Nested functions

% Voltage for all buses in branches array
function [voltFBus, voltTBus] = voltBus_function(branch, bus, ...
    F_BUS, T_BUS, BUS_I, BASE_KV)
% voltFBus=arrayfun(@(x) bus(find(x==bus(:,BUS_I),1),BASE_KV), branch(:,F_BUS));
% voltTBus=arrayfun(@(x) bus(find(x==bus(:,BUS_I),1),BASE_KV), branch(:,T_BUS));    

% voltFBus=zeros(size(branch,1),1);
% voltTBus=zeros(size(branch,1),1);
% for i=1:size(branch,1)
%     voltFBus(i,1)=bus(find(branch(i,F_BUS)==bus(:,BUS_I),1),BASE_KV);
%     voltTBus(i,1)=bus(find(branch(i,T_BUS)==bus(:,BUS_I),1),BASE_KV);
% end

% % % busInBr=ismember(bus(:,BUS_I),branch(:,F_BUS));
% % % voltFBus_aux=bus(busInBr,BASE_KV);
[~,fBusInBr_idx]=ismember(branch(:,F_BUS),bus(:,BUS_I));
voltFBus=bus(fBusInBr_idx,BASE_KV);%Indexing with non-boolean vector
[~,tBusInBr_idx]=ismember(branch(:,T_BUS),bus(:,BUS_I));
voltTBus=bus(tBusInBr_idx,BASE_KV);%Indexing with non-boolean vector
end

% Connectivity for all buses in branches array
function [connBRF, connBRT] = connBR_function(branch, F_BUS, T_BUS)
% connBRF=arrayfun(@(x) sum(branch(:,F_BUS)==x), branch(:,F_BUS));
% connBRT=arrayfun(@(x) sum(branch(:,T_BUS)==x), branch(:,T_BUS));

% connBRF=zeros(size(branch,1),1);
% connBRT=zeros(size(branch,1),1);
% for i=1:size(branch,1)
%     connBRF(i,1)=sum(branch(:,F_BUS)==branch(i,F_BUS));
%     connBRT(i,1)=sum(branch(:,T_BUS)==branch(i,T_BUS));
% end

[cnt,cntIdx] = histc(branch(:,F_BUS), unique(branch(:,F_BUS)));
connBRF = cnt(cntIdx);
[cnt,cntIdx] = histc(branch(:,T_BUS), unique(branch(:,T_BUS)));
connBRT = cnt(cntIdx);

end