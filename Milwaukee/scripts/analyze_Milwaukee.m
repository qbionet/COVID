% Analyze Milwaukee
clear all
close all
clc


% style choices
STYLE.colors.dark = {'#ef476f' '#ffd166' '#06d6a0' '#118ab2' '#adb5bd'};
STYLE.colors.light = {'#f7cad0' '#ffea00' '#b7e4c7' '#a2d2ff'};
STYLE.colors.anal  = {'#0091ad' '#8338ec' '#ff006e'};

STYLE.thin    = 1;
STYLE.medium  = 3;
STYLE.thick   = 5; 
STYLE.fontSize = 12;

STYLE.T_end = 360;
STYLE.XLim_year = [0 366];
STYLE.XTick_months = cumsum([0 31 29 31 30 31 30 31 31 30 31 30 31]);
STYLE.XTickLabel_months = {'Jan' '' '' 'Apr' '' '' 'Jul' '' '' 'Oct' '' '' 'Jan'};

[parentdir,~,~]=fileparts(pwd);

% load statistics when ZIPs are combined
STATS = load([parentdir '/statistics/stats_ZIPcombined.mat']);
ZIP   = STATS.ZIP;
group = STATS.group;
[value index] = max( ZIP.caseEst~=0, [], 2);

% load statistics when ZIPs are not combined
STATS_noComb = load([parentdir '/statistics/stats_ZIPoriginal.mat']);
ZIP_noComb   = STATS_noComb.ZIP;
group_noComb = STATS_noComb.group;
index = find(sum(ZIP_noComb.caseEst),1);

% simulation time horizon
t_start = min(index - 1);
t_end   = 360;
sim_span = t_start:t_end;

% Initial condition 
p0 = 1e-3;             
p0 = [1-p0 p0 0 0];
X0 = ZIP.pop(:,t_start).*p0;
X0_noComb = ZIP_noComb.pop(:,t_start).*p0;


%% Maximum Likelihood Estimation for original ZIPs and corresponding statistics
[psi_noComb error_noComb] = Find_psi(sim_span, X0_noComb, ZIP_noComb, group_noComb, STYLE);
[ZIP_noComb, group_noComb] = SpreadStats(sim_span, ZIP_noComb, group_noComb, psi_noComb);
PlotFit(sim_span, X0_noComb, ZIP_noComb, group_noComb, 'N', STYLE);


%% Maximum Likelihood Estimation for combined ZIPs and corresponding statistics
[psi error] = Find_psi(sim_span, X0, ZIP, group, STYLE);
[ZIP, group] = SpreadStats(sim_span, ZIP, group, psi);
error = PlotFit(sim_span, X0, ZIP, group, 'Y', STYLE);
[mean(error(group.ZIPs{1})) mean(error(group.ZIPs{2})) mean(error(group.ZIPs{3})) mean(error(group.ZIPs{4}))]
[std(error(group.ZIPs{1})) std(error(group.ZIPs{2})) std(error(group.ZIPs{3})) std(error(group.ZIPs{4}))]


%% Where do exposures happen
[ZIP_stats group_stats CS] = ExposureSource(sim_span, X0, ZIP, group, psi, STYLE);
PlotHomophilyExposure(sim_span, ZIP, group, CS, STYLE)


%% Group level differences in trip rate, case rate and vulnerability 
groupDiff_beta(group, STYLE)


%% Violin plot for group level matching
PlotViolin(sim_span, X0, ZIP, group, STYLE)


%% Effects of structural changes (Laplace and isolation)
PlotStructure(sim_span, X0, psi, ZIP, group, STYLE)


%% What would have happened if homophily had not changed compared to before March 1
shift_mean = [];
shift_std  = [];
base_group_posRate = [group.posRate group.posRate*group.pop_mean/sum(group.pop_mean)];
base_inequality = max(base_group_posRate)-min(base_group_posRate);

groupSelect = [1:4];
[difference inequality homo_before homo_after] = EliminateHomophilyShift(ZIP,group,sim_span,X0,groupSelect);
shift_mean = [shift_mean ; difference.mean inequality.mean];
shift_std  = [shift_std ; difference.std inequality.std];

homo_beforeafter = [100*homo_before ; 100*homo_after]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Scripts   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate psi
function [psi_t error_t] = Find_psi(sim_span, X0,  ZIP, group, STYLE)
    
    delta_E = 4;            % mean latency period (4 days)
    delta_I = 3.5;          % mean infectious period (3.5 days)
    
    t_start = sim_span(1);
    t_end   = sim_span(end); 

    nodes = size(X0,1);
    
    X(1:nodes,1:4,t_start) = X0; % initial condition
    
    psi_t = zeros(nodes,t_end);
    for t = t_start:t_end
        N_t    = ZIP.pop(:,t);
        if t == t_end
            N_next = ZIP.pop(:,t);
        else
            N_next = ZIP.pop(:,t+1);
        end
        
        
        % current number of individuals in various disease states for each node
        S_t = round(X(:,1,t));
        E_t = round(X(:,2,t));
        I_t = round(X(:,3,t));
        R_t = N_t - (S_t + E_t + I_t);
                
        % current movement matrix for each node
        W_t = ZIP.W(:,:,t);
        case_t = ZIP.caseEst(:,t);
        
        q_t = I_t./N_t;
        lambda0 = psi_t(:,t-1)./ZIP.area.*W_t'*q_t;
        % calculate psi_t estimate for each node
        for i = 1:nodes
            lambda_hat  = (case_t(i)*N_t(i)./S_t(i) + W_t(i,i)*lambda0(i) - W_t(i,:)*lambda0)/W_t(i,i);
            lambda_hat = min(1, max(0, lambda_hat));
            psi_t(i,t) = lambda_hat*ZIP.area(i)/(W_t(:,i)'*q_t); 
        end
              
        lambda_t = psi_t(:,t)./ZIP.area.*W_t'*q_t;
        error_t(:,t) = [case_t-S_t./N_t.*W_t*lambda_t];
        
        % transitions    
        N_SE = case_t;        
        N_EI = E_t/delta_E;
        N_IR = I_t/delta_I;

        S_next = round(N_next./N_t.*(S_t - N_SE));
        E_next = round(N_next./N_t.*(E_t + N_SE - N_EI));
        I_next = round(N_next./N_t.*(I_t + N_EI - N_IR));
        R_next = N_next - (S_next+E_next+I_next);
        
        X(:,1,t+1) = max([zeros(nodes,1) S_next],[],2);
        X(:,2,t+1) = max([zeros(nodes,1) E_next],[],2);
        X(:,3,t+1) = max([zeros(nodes,1) I_next],[],2);
        X(:,4,t+1) = max([zeros(nodes,1) R_next],[],2);

    end
    psi_t(isnan(psi_t)) = 0; 
end

%% Plot fit results
function error = PlotFit(sim_span, X0, ZIP, group, combined, STYLE)   
    [parentdir,~,~]=fileparts(pwd);

    for k = 1:10
        [casePred_sims(:,:,k)] = SimNet(sim_span, X0, ZIP.W, ZIP.psi, ZIP, group);
    end
    casePred = mean(casePred_sims,3);
    for i = 1:size(ZIP.W,1)
        actual = sum(ZIP.caseEst(i,sim_span),2);
        pred   = sum(casePred(i,sim_span),2);
        meanPop(i) = mean(ZIP.pop(i,sim_span),2);
        error(i) = 100*(pred - actual)/meanPop(i);
    end
    
    
    figure('Position',[0 0 1000 300])
    tlt = tiledlayout(2, 5);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';   
    for g = 1:length(group.name)
        nexttile
        plot(sim_span, squeeze(sum(casePred_sims(group.ZIPs{g},sim_span,:),1)),'o','MarkerFaceColor',STYLE.colors.light{g},'MarkerEdgeColor','none','MarkerSize',4)
        hold on
        plot(sim_span, sum(casePred(group.ZIPs{g},sim_span),1),'Color',STYLE.colors.dark{g},'LineWidth',STYLE.thick)
        plot(sim_span, group.caseEst(g,sim_span),'k','LineWidth',STYLE.thin)
        xticks(STYLE.XTick_months)
        xticklabels(STYLE.XTickLabel_months)
        xlim(STYLE.XLim_year)
        ylabel('case rate (#)')
        title(group.name{g})
        xtickangle(0)
        grid on
    end
      
    nexttile([1 2])
    hold on
    for g = 1:length(group.name)
        index = group.ZIPs{g};
        s = scatter(1e-4*meanPop(group.ZIPs{g}),error(group.ZIPs{g}),50,'filled');

        s.MarkerEdgeColor = 'none';
        s.MarkerFaceColor = STYLE.colors.dark{g};
    end
    % axis square
    grid on
    xlim([0 12])
    ylim([-10 60])
    xlabel('population (10k)')
    ylabel('error (%)')
    
    
    nexttile([1 2])
    hold on
    for g = 1:length(group.name)
        index = group.ZIPs{g};
        s = scatter(100-ZIP.hetero(group.ZIPs{g}),error(group.ZIPs{g}),50,'filled');

        s.MarkerEdgeColor = 'none';
        s.MarkerFaceColor = STYLE.colors.dark{g};
    end
    % axis square
    grid on
    xlim([0 100])
    ylim([-10 60])
    xlabel('homophily (%)')
    ylabel('error (%)')
    
    set(gcf,'renderer','painters');
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize) 

    switch combined
        case 'N'
            print(gcf,[parentdir '/figures/fit_noComb.eps'],'-depsc')
        case 'Y'
            print(gcf,[parentdir '/figures/fit.eps'],'-depsc')
    end

    figure('Position',[0 0 600 250])
    tlt = tiledlayout(3, 2);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';   
    for g = 1:length(group.name)
        nexttile
        plot(sim_span, squeeze(sum(casePred_sims(group.ZIPs{g},sim_span,:),1)),'o','MarkerFaceColor',STYLE.colors.light{g},'MarkerEdgeColor','none','MarkerSize',4)
        hold on
        plot(sim_span, sum(casePred(group.ZIPs{g},sim_span),1),'Color',STYLE.colors.dark{g},'LineWidth',STYLE.thick)
        plot(sim_span, group.caseEst(g,sim_span),'k','LineWidth',STYLE.thin)
        xticks(STYLE.XTick_months)
        xticklabels(STYLE.XTickLabel_months)
        xlim(STYLE.XLim_year)
        yticks([0 2000 4000 6000])
        ylabel('case rate (#)')
        title(group.name{g})
        xtickangle(0)
        grid on
    end
    set(gcf,'renderer','painters');
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize) 
    switch combined
        case 'N'
            print(gcf,[parentdir '/figures/fit_compact_noComb.eps'],'-depsc')
        case 'Y'
            print(gcf,[parentdir '/figures/fit_compact.eps'],'-depsc')
    end
end




%% Simulate network
function [casePred fractions caseSource risk X] = SimNet(sim_span, X0, W, psi_t, ZIP, group)
    
    delta_E = 4;            % mean latency period (4 days)
    delta_I = 3.5;          % mean infectious period (3.5 days)
    
    t_start = sim_span(1);
    t_end   = sim_span(end); 

    nodes = size(X0,1);
    % initial condition
    X(1:nodes,1:4,t_start) = X0;
    case_all   = zeros(nodes,t_end);
    case_self  = zeros(nodes,t_end);
    case_group = zeros(4,t_end);


    for t = t_start:t_end
        N_t    = ZIP.pop(:,t);
        if t == t_end
            N_next = ZIP.pop(:,t);
        else
            N_next = ZIP.pop(:,t+1);
        end
               
        % current number of individuals in various disease states for each node
        S_t = round(X(:,1,t));
        E_t = round(X(:,2,t));
        I_t = round(X(:,3,t));
        R_t = N_t - (S_t + E_t + I_t);
                
        % current movement matrix for each node
        W_t = squeeze(W(:,:,t));
               
        q_t = I_t./N_t;
        
        for i = 1:length(q_t)
            Wq_sum = 0;
            for j = 1:length(q_t)
               Wq_sum = Wq_sum + W_t(j,i)*q_t(j); 
            end
            lambda_t(i,1) = psi_t(i,t)/ZIP.area(i)*Wq_sum;          
        end
        lambda_t(lambda_t>1) = 1;
        

        % transitions
        N_SE = poissrnd(S_t./N_t.*W_t*lambda_t);       
        N_EI = poissrnd(E_t/delta_E);
        N_IR = poissrnd(I_t/delta_I);

        S_next = round(N_next./N_t.*(S_t - N_SE));
        E_next = round(N_next./N_t.*(E_t + N_SE - N_EI));
        I_next = round(N_next./N_t.*(I_t + N_EI - N_IR));
        R_next = N_next - (S_next+E_next+I_next);
        
        X(:,1,t+1) = max([zeros(nodes,1) S_next],[],2);
        X(:,2,t+1) = max([zeros(nodes,1) E_next],[],2);
        X(:,3,t+1) = max([zeros(nodes,1) I_next],[],2);
        X(:,4,t+1) = max([zeros(nodes,1) R_next],[],2);

        casePred(:,t) = N_SE;
        L(:,t) = lambda_t;
        Q(:,t) = q_t;
        
        
        %% Where do infections come from?
        % n x n, each row contains the % of exposures from all nodes
        fracMatrix = W_t.*repmat(lambda_t,1,size(W_t,2));
        fracMatrix = fracMatrix./repmat(sum(fracMatrix,2),1,size(fracMatrix,2));
        
        % n x n, each row contains the # of exposures from all nodes
        caseSource(:,:,t) = repmat(N_SE,1,size(fracMatrix,2)).*fracMatrix;
                
        case_all(:,t) = N_SE;
        case_all(isnan(case_all)) = 0;
        for g = 1:length(group.name)
            index_group = group.ZIPs{g};
            for i = 1:length(group.ZIPs{g})               
                index_self = group.ZIPs{g}(i);
                rate_all   = W_t(index_self,:)*lambda_t;
                rate_group = W_t(index_self,index_group)*lambda_t(index_group);
                rate_self  = W_t(index_self,index_self)*lambda_t(index_self);
                                
                case_self(index_self,t)  = rate_self/rate_all*case_all(index_self,t);
                case_group(index_self,t) = rate_group/rate_all*case_all(index_self,t);
            end
        end
        
        for g1 = 1:length(group.name)
            index1 = group.ZIPs{g1};
            for g2 = 1:length(group.name)
                index2 = group.ZIPs{g2};
                risk.group(g1,g2,t) = sum(W_t(index1,index2)*lambda_t(index2))/sum(W_t(index1,:)*lambda_t(:));              
            end
        end
        for i = 1:nodes
            for g = 1:length(group.name)
                index = group.ZIPs{g};
                risk.selfGroup(i,g,t) = sum(W_t(i,index)*lambda_t(index))/sum(W_t(i,:)*lambda_t(:));
            end
        end
        
        for i = 1:nodes
           risk.self(i,t) = W_t(i,i)*lambda_t(i)/(W_t(i,:)*lambda_t(:)); 
        end
        
    end
 
    % fraction from own node
    fractions.ZIP.self  = 100*sum(case_self,2,'omitnan')./sum(case_all,2,'omitnan');
    
    % fraction from own demographic
    fractions.ZIP.demo = 100*sum(case_group,2,'omitnan')./sum(case_all,2,'omitnan');
    for g = 1:length(group.name)
        index_group = group.ZIPs{g};
        % fraction from own node
        fractions.group.self(g,1) = 100*sum(case_self(index_group,:),[1 2],'omitnan')/sum(case_all(index_group,:),[1 2],'omitnan');
        
        % fraction from own demographic
        fractions.group.demo(g,1) = 100*sum(case_group(index_group,:),[1 2],'omitnan')/sum(case_all(index_group,:),[1 2],'omitnan');
    end
end

%% Compute additional group and ZIP stats
function [ZIP, group] = SpreadStats(sim_span,ZIP,group,psi_t)

    for i = 1:length(ZIP.name)
        ZIP.popDens(i,sim_span)  = ZIP.pop(i,sim_span)/ZIP.area(i);
        ZIP.selfRate(i,sim_span) = squeeze(ZIP.W(i,i,sim_span))'./ZIP.pop(i,sim_span);

        ZIP.beta(i,sim_span) = psi_t(i,sim_span).*ZIP.popDens(i,sim_span).*ZIP.selfRate(i,sim_span).^2;
        ZIP.psi(i,sim_span)  = psi_t(i,sim_span);
        ZIP.posRate(i) = 100*sum(ZIP.caseEst(i,sim_span),2)/mean(ZIP.pop(i,sim_span),2);
        ZIP.hetero(i) = 100 - 100*sum(ZIP.W(i,i,sim_span),3)/sum(ZIP.W(i,:,sim_span),[2 3]);
        ZIP.beta_mean(i) = mean(ZIP.beta(i,sim_span),2);
        ZIP.pop_mean(i) = mean(ZIP.pop(i,sim_span),2);
        ZIP.vuln(i,sim_span) = ZIP.psi(i,sim_span).*ZIP.popDens(i,sim_span);
        ZIP.vuln_mean(i) = mean(ZIP.vuln(i,sim_span),2);
        ZIP.tripRate_mean(i) = mean(ZIP.tripRate(i,sim_span),2);
        
        
        ZIP.beta1(i,sim_span) = ZIP.psi(i,sim_span).*ZIP.popDens(i,sim_span);
        ZIP.beta2(i,sim_span) = ZIP.tripRate(i,sim_span).^2;
        ZIP.betaRatio(i,sim_span) = ZIP.beta2(i,sim_span)./ZIP.beta1(i,sim_span);
        
        ZIP.beta1_mean(i) = mean(ZIP.beta1(i,sim_span),2);
        ZIP.beta2_mean(i) = mean(ZIP.beta2(i,sim_span),2);
        ZIP.betaRatio_mean(i) = mean(ZIP.betaRatio(i,sim_span),2,'omitnan');
        
    end

    for g = 1:length(group.name)
        index = group.ZIPs{g};
        group.popDens(g,sim_span)  = group.pop(g,sim_span)./sum(ZIP.area(index));
        group.selfRate(g,sim_span) = squeeze(sum(ZIP.W(index,index,sim_span),[1,2]))'./group.pop(g,sim_span);
        group.psi(g,sim_span) = sum(ZIP.psi(index,sim_span).*ZIP.pop(index,sim_span),1)./sum(ZIP.pop(index,sim_span),1);
        group.beta(g,sim_span) = sum(ZIP.beta(index,sim_span).*ZIP.pop(index,sim_span),1)./sum(ZIP.pop(index,sim_span),1);

        
        group.vuln(g,sim_span) = sum(ZIP.vuln(index,sim_span).*ZIP.pop(index,sim_span),1)./sum(ZIP.pop(index,sim_span),1);
        
        group.posRate(g) = 100*sum(group.caseEst(g,sim_span),2)/mean(group.pop(g,sim_span),2);
        
        group.beta1(g,sim_span) = sum(ZIP.beta1(index,sim_span).*ZIP.pop(index,sim_span),1)./sum(ZIP.pop(index,sim_span),1);
        group.beta2(g,sim_span) = sum(ZIP.beta2(index,sim_span).*ZIP.pop(index,sim_span),1)./sum(ZIP.pop(index,sim_span),1);
        group.betaRatio(g,sim_span) = sum(ZIP.betaRatio(index,sim_span).*ZIP.pop(index,sim_span),1)./sum(ZIP.pop(index,sim_span),1);
        
        group.hetero(g) = 100 - 100*sum(ZIP.W(index,index,sim_span),[1 2 3])/sum(ZIP.W(index,:,sim_span),[1 2 3]);
        group.hetero_t(g,sim_span) = 100 - 100*sum(ZIP.W(index,index,sim_span),[1 2])./sum(ZIP.W(index,:,sim_span),[1 2]);
    end
    group.tripRate_mean = mean(group.tripRate(:,sim_span),2);
    group.popDens_mean  = mean(group.popDens(:,sim_span),2);
    group.pop_mean      = mean(group.pop(:,sim_span),2);
    group.psi_mean      = mean(group.psi(:,sim_span),2);
    group.beta_mean     = mean(group.beta(:,sim_span),2);
    group.vuln_mean        = mean(group.vuln(:,sim_span),2);
    
    group.beta1_mean = mean(group.beta1(:,sim_span),2);
    group.beta2_mean = mean(group.beta2(:,sim_span),2);
    group.betaRatio_mean = mean(group.betaRatio(:,sim_span),2);
end



%% Trip risk
function [ZIP_stats group_stats CS] = ExposureSource(sim_span, X0, ZIP, group, psi_t, STYLE)
    [parentdir,~,~]=fileparts(pwd);

    K_rep = 100; % number of replicates

    figure('Position',[0 0 1000 250])
    tlt = tiledlayout(1, 5);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';  


    % Compute where exposures happen
    for k = 1:K_rep
        [casePred fractions(k) cs risk] = SimNet(sim_span, X0, ZIP.W, psi_t, ZIP, group);


        % cs is casesource: n x n x t, each row contains the # of exposures
        % from all nodes at time t

        CS{k} = cs;
        selfRisk(:,:,k) = risk.self;
        for g = 1:length(group.name)
            groupRisk(g,:,k) = squeeze(risk.group(g,g,:));
        end
        
        for g = 1:length(group.name)
            index = group.ZIPs{g};
            selfGroupRisk(g,:,k) = squeeze(mean(risk.selfGroup(index,g,:),1));
        end
        
    end
    
    sourceNode_mean = zeros(length(fractions(1).ZIP.self),1);
    sourceGroup_mean = zeros(length(fractions(1).ZIP.self),1);
    for k = 1:K_rep
        sourceNode_mean  = sourceNode_mean + fractions(k).ZIP.self;
        sourceGroup_mean = sourceGroup_mean + fractions(k).ZIP.demo;
    end
    sourceNode_mean = sourceNode_mean/K_rep;
    sourceGroup_mean = sourceGroup_mean/K_rep;
    
    selfRisk_mean  = mean(selfRisk,3,'omitnan');
    groupRisk_mean = mean(groupRisk,3,'omitnan');
    selfGroupRisk_mean = mean(selfGroupRisk,3,'omitnan');

    % Compute self trips
    for i = 1:size(ZIP.W,1)
        w_self  = squeeze(ZIP.W(i,i,:));
        w_total = squeeze(sum(ZIP.W(i,:,:),2));
        selfTrips(i,:) = w_self./w_total;
    end
    
    % Compute self group trips
    for g = 1:length(group.name)
        index = group.ZIPs{g};       
        w_self  = squeeze(sum(ZIP.W(index,index,:),[1 2]));
        w_total = squeeze(sum(ZIP.W(index,:,:),[1 2]));
        groupTrips(g,:) = w_self./w_total;
    end
    
    % Self trips
    nexttile
    hold on
    for g = 1:length(group.name)
        index = group.ZIPs{g};
        selfTrips_groupMean(g,:) = mean(selfTrips(index,:),1);
    end
    for g = 1:length(group.name)
        plot(100*selfTrips_groupMean(g,:),'LineWidth',STYLE.medium,'Color',STYLE.colors.dark{g})
    end
    axis square
    grid on
    ylim([0 100])
    xlim(STYLE.XLim_year)
    xticks(STYLE.XTick_months)
    xticklabels(STYLE.XTickLabel_months)
    xtickangle(0)
    ylabel('homophily (%)')
    title('trips within own ZIP')
    
    
    % Self risk
    nexttile
    hold on
    for g = 1:length(group.name)
        index = group.ZIPs{g};
        selfRisk_groupMean(g,:) = mean(selfRisk_mean(index,:),1);
    end
    selfRisk_groupMean(find(selfRisk_groupMean == 0)) = NaN;
    for g = 1:4
        plot(100*selfRisk_groupMean(g,:),'LineWidth',STYLE.medium,'Color',STYLE.colors.dark{g})
    end
    axis square
    grid on
    ylim([0 100])
    xlim(STYLE.XLim_year)
    xticks(STYLE.XTick_months)
    xticklabels(STYLE.XTickLabel_months)
    xtickangle(0)
    ylabel('expected probability (%)')
    title('exposures within own ZIP')
    
    % Group trips
    nexttile
    hold on
    for g = 1:length(group.name)
        plot(100*groupTrips(g,:),'LineWidth',STYLE.medium,'Color',STYLE.colors.dark{g})
    end
    axis square
    grid on
    ylim([0 100])
    xlim(STYLE.XLim_year)
    xticks(STYLE.XTick_months)
    xticklabels(STYLE.XTickLabel_months)
    xtickangle(0)
    ylabel('homophily (%)')
    title('trips within own group')
    
    % group risk
    nexttile
    hold on
    
    groupRisk_mean(find(groupRisk_mean == 0)) = NaN;
    selfGroupRisk_mean(find(selfGroupRisk_mean == 0)) = NaN;
    for g = 1:length(group.name)
        plot(100*selfGroupRisk_mean(g,:),'LineWidth',STYLE.medium,'Color',STYLE.colors.dark{g})
    end
    axis square
    grid on
    ylim([0 100])
    xlim(STYLE.XLim_year)
    xticks(STYLE.XTick_months)
    xticklabels(STYLE.XTickLabel_months)
    xtickangle(0)
    ylabel('expected probability (%)')
    title('exposures within own group')
    
    
    stats = [mean(selfTrips_groupMean,2,'omitnan')  mean(selfTrips_groupMean,2,'omitnan')  mean(groupTrips,2,'omitnan')  mean(selfGroupRisk_mean,2,'omitnan')];
    [row,column] = size(stats)
    for i = 1:column
        stats(row+1,i) = group.pop_mean'*stats(1:row,i)/sum(group.pop_mean);
    end
    
    
    for g = 1:length(group.name)
        for i = 1:length(group.ZIPs{g})
            group_index = group.ZIPs{g};
            ZIP_index = group_index(i);
                        
            ZIP_stats{g}.name(i)         = ZIP.name(ZIP_index);
            ZIP_stats{g}.node_hetero(i)  = 100*sum(ZIP.W(ZIP_index,ZIP_index,:),[1:3])/sum(ZIP.W(ZIP_index,:,:),[1:3]);
            ZIP_stats{g}.group_hetero(i) = 100*sum(ZIP.W(ZIP_index,group_index,:),[1:3])/sum(ZIP.W(ZIP_index,:,:),[1:3]);
            
            ZIP_stats{g}.node_source(i)  = sourceNode_mean(ZIP_index);
            ZIP_stats{g}.group_source(i) = sourceGroup_mean(ZIP_index);
                      
            ZIP_stats{g}.summary(i,:) = [ZIP_stats{g}.name(i) ZIP.pop_mean(ZIP_index) ZIP.posRate(ZIP_index) ZIP_stats{g}.node_hetero(i) ZIP_stats{g}.group_hetero(i) ZIP_stats{g}.node_source(i) ZIP_stats{g}.group_source(i)];
        end

        group_stats{g}.ZIP_exposure    = mean(100*selfRisk_groupMean(g,sim_span),'omitnan');
        group_stats{g}.ZIP_homophily   = mean(100*selfTrips_groupMean(g,sim_span),'omitnan');

        group_stats{g}.group_exposure  = mean(100*selfGroupRisk_mean(g,sim_span),'omitnan');
        group_stats{g}.group_homophily = mean(100*groupTrips(g,sim_span),'omitnan');
    end
    
 
    set(gcf,'renderer','painters');
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize) 
    print(gcf,[parentdir '/figures/exposureSource.eps'],'-depsc')
    
end



%% Radial plots for homophily and sources of exposure
function PlotHomophilyExposure(sim_span, ZIP, group, CS, STYLE)
    [parentdir,~,~]=fileparts(pwd);

    circr = @(radius,rad_ang)  [radius*cos(rad_ang);  radius*sin(rad_ang)];
    N_rad = 1e4;
    
    
    time_span = sim_span(1):sim_span(end);
    % average ZIP homophily
    total_ZIP_self_trip_city = 0;
    total_group_self_trip_city = 0;
    for g = 1:length(group.name)
        index = group.ZIPs{g};
        total_group_trip    = sum(ZIP.W(index,:,time_span),[1 2 3],'omitnan');
        total_ZIP_self_trip = 0;
        for i = 1:length(index);
            total_ZIP_self_trip = total_ZIP_self_trip + sum(squeeze(ZIP.W(index(i),index(i),time_span)),'omitnan');
        end
        ZIP_homophily(g) = 100*total_ZIP_self_trip/total_group_trip;
        total_ZIP_self_trip_city = total_ZIP_self_trip_city + total_ZIP_self_trip;
    
        for i = 1:length(group.name)
            target = group.ZIPs{i};
            total_trip_to_target_group = sum(ZIP.W(index,target,time_span),[1 2 3],'omitnan');
            group_homophily(g,i) = 100*total_trip_to_target_group/total_group_trip;
        end
        total_group_self_trip = sum(ZIP.W(index,index,time_span),[1 2 3],'omitnan');
        total_group_self_trip_city = total_group_self_trip_city + total_group_self_trip;
        
    end
    
    total_city_trip = sum(ZIP.W(:,:,time_span),[1 2 3],'omitnan');
    city_homophily = 100*[total_ZIP_self_trip_city total_group_self_trip_city]/total_city_trip;
    
    
    for k = 1:length(CS)
        cs = CS{k};
        % average ZIP homophily
        total_ZIP_self_case_city(k) = 0;
        total_group_self_case_city(k) = 0;
        for g = 1:length(group.name)
            index = group.ZIPs{g};
            total_group_case    = sum(cs(index,:,sim_span),[1 2 3],'omitnan');
            total_ZIP_self_case = 0;
            for i = 1:length(index);
                total_ZIP_self_case = total_ZIP_self_case + sum(squeeze(cs(index(i),index(i),sim_span)),'omitnan');
            end
            ZIP_case_rep(k,g) = 100*total_ZIP_self_case/total_group_case;
            total_ZIP_self_case_city(k) = total_ZIP_self_case_city(k) + total_ZIP_self_case;
    
            for i = 1:length(group.name)
                target = group.ZIPs{i};
                total_case_from_target_group = sum(cs(index,target,sim_span),[1 2 3],'omitnan');
                group_case_rep(k,g,i) = 100*total_case_from_target_group/total_group_case;
            end
            total_group_self_case = sum(cs(index,index,sim_span),[1 2 3],'omitnan');
            total_group_self_case_city(k) = total_group_self_case_city(k) + total_group_self_case;
        end
        total_city_case = sum(cs(:,:,sim_span),[1 2 3],'omitnan');
    end
    ZIP_case   = mean(ZIP_case_rep,1,'omitnan');
    group_case = squeeze(mean(group_case_rep,1,'omitnan'));
    city_homophily = 100*[mean(total_ZIP_self_case_city./total_city_case) mean(total_group_self_trip_city./total_city_trip)];
    
    
    figure('Position',[0 0 800 400])
    tlt = tiledlayout(2, 5);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'none';
    
    for i = 1:length(group.name)
        nexttile(i)
        hold on
    
        N = 20;
        for k = 0:N-1
            rho = 1.4;
            theta = pi/2-k*2*pi/N;
            [x y] = pol2cart(theta,rho);
            plot([0 x],[0 y],'k','LineWidth',0.2)
        end
    
        N = 5;
        for k = 0:N-1
            rho = 1.7;
            theta = pi/2-k*2*pi/N;
            [x y] = pol2cart(theta,rho);
            text(x,y,[num2str(100/N*k) '%'],'HorizontalAlignment','center')
        end
        
        data = cumsum(circshift(group_case(i,:),-i+1));
        color_code = circshift([1:4],-i+1);
        
        for g = length(group.name):-1:1
            radius = 1; 
            r_angl = pi/2-linspace(0, data(g), N_rad)/100*2*pi;  
            xy_r = circr(radius,r_angl);
            plot(xy_r(1,:), xy_r(2,:),'LineWidth',20,'Color',hex2rgb(STYLE.colors.dark{color_code(g)}))
        end
    
        radius = 0.45; 
        r_angl = sum(group_case(i,1:i-1))/100*2*pi+pi/2-linspace(sum(group_case(i,1:i-1)), sum(group_case(i,1:i-1))+ZIP_case(i), N_rad)/100*2*pi;  
        xy_r = circr(radius,r_angl);
        plot(xy_r(1,:), xy_r(2,:),'LineWidth',20,'Color',hex2rgb(STYLE.colors.light{i}))
        axis square
        xlim([-2 2]), ylim([-2 2])
        xticks([]), yticks([])
        title(group.name{i})
        set(gca,'Visible','off')
    end
    
    for i = 1:length(group.name)
        nexttile(5+i)
        hold on
    
        N = 20;
        for k = 0:N-1
            rho = 1.4;
            theta = pi/2-k*2*pi/N;
            [x y] = pol2cart(theta,rho);
            plot([0 x],[0 y],'k','LineWidth',0.2)
        end
    
        N = 5;
        for k = 0:N-1
            rho = 1.7;
            theta = pi/2-k*2*pi/N;
            [x y] = pol2cart(theta,rho);
            text(x,y,[num2str(100/N*k) '%'],'HorizontalAlignment','center')
        end
    
        data = cumsum(circshift(group_homophily(i,:),-i+1));
        color_code = circshift([1:length(group.name)],-i+1);
        
        for g = length(group.name):-1:1
            radius = 1; 
            r_angl = pi/2-linspace(0, data(g), N_rad)/100*2*pi;  
            xy_r = circr(radius,r_angl);
            plot(xy_r(1,:), xy_r(2,:),'LineWidth',20,'Color',hex2rgb(STYLE.colors.dark{color_code(g)}))
        end
    
        radius = 0.45; 
        r_angl = sum(group_homophily(i,1:i-1))/100*2*pi+pi/2-linspace(sum(group_homophily(i,1:i-1)), sum(group_homophily(i,1:i-1))+ZIP_homophily(i), N_rad)/100*2*pi;  
        xy_r = circr(radius,r_angl);
        plot(xy_r(1,:), xy_r(2,:),'LineWidth',20,'Color',hex2rgb(STYLE.colors.light{i}))
        axis square
        xlim([-2 2]), ylim([-2 2])
        xticks([]), yticks([])
        title(group.name{i})
        set(gca,'Visible','off')
    
    end
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize) 
    
    set(gcf,'renderer','opengl');
    print(gcf,[parentdir '/figures/homophily_exposure_circles.eps'],'-depsc','-r1200')

end


%% Group differences
function groupDiff_beta(group, STYLE)
 
    [parentdir,~,~]=fileparts(pwd);
    
  
    figure('Position',[0 0 250 350])
    tlt = tiledlayout(3, 1);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';

    nexttile
    b = barh(group.posRate/min(group.posRate),'FaceColor','flat','EdgeColor','none');
    for g = 1:length(group.name)
        b.CData(g,:) = hex2rgb(STYLE.colors.dark{g});
    end
    grid on
    ylim([0.5 length(group.name)+0.5])
    xlim([0 3])
    xticks([0:0.5:3])
    xticklabels({0 '' 1 '' 2 '' 3})
    yticklabels({'Black' 'Latinx' 'White' 'Mixed'})
    title('case rate')
    xtickangle(0)

    nexttile
    b = barh(group.beta1_mean/min(group.beta1_mean),'FaceColor','flat','EdgeColor','none');
    for g = 1:length(group.name)
        b.CData(g,:) = hex2rgb(STYLE.colors.dark{g});
    end
    grid on
    ylim([0.5 length(group.name)+0.5])
    xlim([0 3])
    xticks([0:0.5:3])
    xticklabels({0 '' 1 '' 2 '' 3})
    yticklabels({'Black' 'Latinx' 'White' 'Mixed'})
    title('vulnerability')
    xtickangle(0)
    
    nexttile
    b = barh(group.tripRate_mean/min(group.tripRate_mean),'FaceColor','flat','EdgeColor','none');
    for g = 1:length(group.name)
        b.CData(g,:) = hex2rgb(STYLE.colors.dark{g});
    end
    grid on
    ylim([0.5 length(group.name)+0.5])
    xlim([0 3])
    xticks([0:0.5:3])
    xticklabels({0 '' 1 '' 2 '' 3})
    yticklabels({'Black' 'Latinx' 'White' 'Mixed'})
    title('trip rate')
    xtickangle(0)

    set(gcf,'renderer','painters');
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize) 
    print(gcf,[parentdir '/figures/groupDiff.eps'],'-depsc')
end





%% Calculate hypothetical statistics of ZIPS and groups
function hypoStats = HypoStats(sim_span, X0, W, psi, ZIP, group)

    N_rep = 10;
    
    casePred = [];
    for k = 1:N_rep
        casePred(:,:,k) = SimNet(sim_span, X0, W, psi, ZIP, group);
    end


    % calculate case count
    casePred_mean = mean(casePred,3);   % mean exposure across runs
    hypoStats.ZIP.posRate = 100*sum(casePred_mean,2,'omitnan')./mean(ZIP.pop,2,'omitnan');
    for g = 1:length(group.name)
        group_pos(g,:)  = sum(casePred_mean(group.ZIPs{g},:),1,'omitnan');
    end
    hypoStats.group.posRate  = 100*sum(group_pos,2)./mean(group.pop,2);
    hypoStats.group.posCum   = 100*cumsum(group_pos,2)./mean(group.pop,2);

    % calculate trip rate for ZIPs and groups
    ZIP_trip_hypo = squeeze(sum(W,2));
    hypoStats.ZIP.tripRate = ZIP_trip_hypo./ZIP.pop;
    dW = squeeze(sum(W,2));
    for g = 1:length(group.name)
        group_trip_hypo(g,:) = sum(dW(group.ZIPs{g},:),1);
    end
    hypoStats.group.tripRate = group_trip_hypo./group.pop;

    
    % calculate pop dens for ZIPS and groups
    for i = 1:length(ZIP.area)
        hypoStats.ZIP.popDens(i,:) = ZIP.pop(i,:)./ZIP.area(i);
    end
    for g = 1:length(group.name)
        index = group.ZIPs{g};
        hypoStats.group.popDens(g,:)  = group.pop(g,:)./sum(ZIP.area(index));
    end
    
    
    % calculate psi for ZIPS and groups
    hypoStats.ZIP.psi = psi;
    for g = 1:length(group.name)
        index = group.ZIPs{g};
        hypoStats.group.psi(g,:) = sum(psi(index,:).*ZIP.pop(index,:),1)./sum(ZIP.pop(index,:),1);
    end
    
    % calculate vulnerability for ZIPS and groups
    hypoStats.ZIP.vuln = psi.*hypoStats.ZIP.popDens;
    for g = 1:length(group.name)
        index = group.ZIPs{g};
        hypoStats.group.vuln(g,:) = sum(hypoStats.ZIP.vuln(index,:).*ZIP.pop(index,:),1)./sum(ZIP.pop(index,:),1);
    end
    
    
    % calculate beta for ZIPs and groups
    for i = 1:length(ZIP.area)
        ZIP_selfRate(i,:) = squeeze(W(i,i,:))'./ZIP.pop(i,:);
        hypoStats.ZIP.beta(i,:) = hypoStats.ZIP.psi(i,:).*hypoStats.ZIP.popDens(i,:).*ZIP_selfRate(i,:).^2;
    end
    for g = 1:length(group.name)
        index = group.ZIPs{g};
        group_selfRate(g,:) = squeeze(sum(W(index,index,:),[1,2]))'./group.pop(g,:);
        hypoStats.group.beta(g,:) = sum(hypoStats.ZIP.beta(index,:).*ZIP.pop(index,:),1)./sum(ZIP.pop(index,:),1);
    end
    
end



%% Matching scatter plots
function rateDiff = MatchingScatters(sim_span,X0,ZIP,group,match,STYLE)
    
    rateDiff(1) = group.posRate(match.target) - group.posRate(match.source);

    ZIPs_toRescale = group.ZIPs{match.source};

    %% rescale triprate mean and temporal
    W_hypo    = ZIP.W;
    psi_hypo  = ZIP.psi;   
    for i = 1:length(ZIPs_toRescale)
        index = ZIPs_toRescale(i);
        factor_t = group.tripRate(match.target,:)./ZIP.tripRate(index,:);
        factor_t(isnan(factor_t)) = 0;
        tripRate_old = [squeeze(sum(W_hypo(index,:,:),2))]'./ZIP.pop(index,:);
        for t = 1:length(factor_t)
            W_hypo(index,:,t) = factor_t(t)*W_hypo(index,:,t);
        end
        tripRate_new = [squeeze(sum(W_hypo(index,:,:),2))]'./ZIP.pop(index,:);
        factor_mean = mean(tripRate_old(sim_span),'omitnan')/mean(tripRate_new(sim_span),'omitnan');
        W_hypo(index,:,:) = factor_mean*W_hypo(index,:,:);
    end
    factor_W = group.tripRate_mean(match.target)/group.tripRate_mean(match.source);
    W_hypo(ZIPs_toRescale,:,:) = factor_W*W_hypo(ZIPs_toRescale,:,:);
    hypoStats = HypoStats(sim_span, X0, W_hypo, psi_hypo, ZIP, group); 
    rateDiff(2) = hypoStats.group.posRate(match.target) - hypoStats.group.posRate(match.source);
    
    %% rescale psi mean and temporal
    W_hypo    = ZIP.W;
    psi_hypo  = ZIP.psi;    
    for i = 1:length(ZIPs_toRescale)
        index = ZIPs_toRescale(i);
        factor_t = group.vuln(match.target,:)./ZIP.vuln(index,:);
        psi_old = psi_hypo(index,:);
        vuln_old = ZIP.vuln(index,:);
        for t = 1:length(factor_t)
            psi_hypo(index,t) = factor_t(t)*psi_hypo(index,t);
            vuln_new(t) = factor_t(t)*vuln_old(t);
        end
        psi_new = psi_hypo(index,:);
        factor_mean = mean(vuln_old(sim_span),'omitnan')/mean(vuln_new(sim_span),'omitnan');
        psi_hypo(index,:) = factor_mean*psi_hypo(index,:);
    end
    psi_hypo(isnan(psi_hypo)) = 0;        
    factor_psi  = group.vuln_mean(match.target)/group.vuln_mean(match.source);
    psi_hypo(ZIPs_toRescale,:) = factor_psi*psi_hypo(ZIPs_toRescale,:);    
    hypoStats = HypoStats(sim_span, X0, W_hypo, psi_hypo, ZIP, group);
    rateDiff(3) = hypoStats.group.posRate(match.target) - hypoStats.group.posRate(match.source);
 
    % beta matching (mean matching via trip rate)    
    W_hypo    = ZIP.W;
    psi_hypo  = ZIP.psi;

    % rescale tripRate temporal
    for i = 1:length(ZIPs_toRescale)
        index = ZIPs_toRescale(i);
        factor_t = group.tripRate(match.target,:)./ZIP.tripRate(index,:);
        factor_t(isnan(factor_t)) = 0;
        tripRate_old = [squeeze(sum(W_hypo(index,:,:),2))]'./ZIP.pop(index,:);
        for t = 1:length(factor_t)
            W_hypo(index,:,t) = factor_t(t)*W_hypo(index,:,t);
        end
        tripRate_new = [squeeze(sum(W_hypo(index,:,:),2))]'./ZIP.pop(index,:);
        factor_mean = mean(tripRate_old(sim_span))/mean(tripRate_new(sim_span));
        W_hypo(index,:,:) = factor_mean*W_hypo(index,:,:);
    end
    
    % rescale psi temporal
    for i = 1:length(ZIPs_toRescale)
        index = ZIPs_toRescale(i);
        factor_t = group.vuln(match.target,:)./ZIP.vuln(index,:);
        psi_old = psi_hypo(index,:);
        vuln_old = ZIP.vuln(index,:);
        for t = 1:length(factor_t)
            psi_hypo(index,t) = factor_t(t)*psi_hypo(index,t);
            vuln_new(t) = factor_t(t)*vuln_old(t);
        end
        psi_new = psi_hypo(index,:);
        factor_mean = mean(vuln_old(sim_span),'omitnan')/mean(vuln_new(sim_span),'omitnan');
        psi_hypo(index,:) = factor_mean*psi_hypo(index,:);
        psi_hypo(index,isnan(psi_hypo(index,:))) = 0;
    end;
    psi_hypo(isnan(psi_hypo)) = 0;

    % mean of target beta
    beta_target = group.beta_mean(match.target);

    % modified beta (slight change due to trajectory shifts
    for i = 1:length(ZIPs_toRescale)
        index = ZIPs_toRescale(i);

        ZIP_popDens(index,:)  = ZIP.pop(index,:)/ZIP.area(index);
        ZIP_selfRate(index,:) = squeeze(W_hypo(index,index,:))'./ZIP.pop(index,:);

        ZIP_beta(index,:) = psi_hypo(index,:).*ZIP_popDens(index,:).*ZIP_selfRate(index,:).^2;
    end
    group_beta(:) = sum(ZIP_beta(ZIPs_toRescale,:).*ZIP.pop(ZIPs_toRescale,:),1)./sum(ZIP.pop(ZIPs_toRescale,:),1);
    beta_source = mean(group_beta(sim_span),2);

    matchType = 2;
    % how to match mean beta
    switch matchType
        case 0
            factor_psi = 1;
            factor_W   = 1;
            TITLE = 'trajectory matching';
        case 1 
            factor_psi = beta_target/beta_source;
            factor_W   = 1;
            TITLE = 'mean matching via vulnerability';
        case 2
            factor_psi = 1;
            factor_W   = sqrt(beta_target/beta_source);
            TITLE = 'mean matching via trip rate';
    end
    psi_hypo(ZIPs_toRescale,:) = factor_psi*psi_hypo(ZIPs_toRescale,:);
    W_hypo(ZIPs_toRescale,:,:) = factor_W*W_hypo(ZIPs_toRescale,:,:);

    hypoStats = HypoStats(sim_span, X0, W_hypo, psi_hypo, ZIP, group);
    rateDiff(4) = hypoStats.group.posRate(match.target) - hypoStats.group.posRate(match.source);
    
end

%% Violin plot for matchings
function PlotViolin(sim_span, X0, ZIP, group, STYLE)
    [parentdir,~,~]=fileparts(pwd);
    
    for i = 1:length(group.name)
        for j = setdiff([1:length(group.name)],i)
        match.source = i;
        match.target = j;
        rateDiff = MatchingScatters(sim_span,X0,ZIP,group,match,STYLE);
        RATEDIFF(i,j,:) = rateDiff;
        end
    end

    for k = 1:4
        H{k} = squeeze(abs(RATEDIFF(:,:,k)))
    end
    
    data  = [];
    count = 0;
    scenarios = {'baseline' ; 'trip rate'; 'vulnerability'; 'approx tx rate'}; 
    ViolinColor = 0.2*ones(4,3);
    for i = 1:4
        h = unique(H{i}(find(H{i}>0)));
        S(i,:) = [min(h) max(h) mean(h) std(h)];
    
        data = [data ; h];
        for j = 1:length(h)
            count = count+1;
            name{count,1} = scenarios{i};
        end
        
    end
    figure('Position',[0 0 600 250])
    violinplot(data, cellstr(name),'ViolinColor',ViolinColor,'GroupOrder',scenarios,'BoxColor',[0 0 0]);
    ylabel('case rate difference between groups (%.)');
    xlim([0.5 4.5]), % ylim([0, 30]), yticks([0:10:30])
    set(gca,'YGrid','on')
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize) 
        
    set(gcf,'renderer','opengl');
    print(gcf,[parentdir '/figures/matching_violin.eps'],'-depsc','-r1200')
end













%% Laplace smoothing
function [W_L heterophily] = LaplaceSmoothing_group(W_base,alpha,group_index,group)
    N_nodes = size(W_base,2);
    
    for t = 1:size(W_base,3)
        W_L(:,:,t) = W_base(:,:,t);
        index = group.ZIPs{group_index};
        w = W_base(index,:,t);
        d = repmat(sum(w,2),1,N_nodes);
        w_new = (w+alpha*d.^2/N_nodes)./(1+alpha*d);
        W_L(index,:,t) = w_new;
    end
    
    for g = 1:length(group.name)
        index = group.ZIPs{g};
        heterophily(g) = 100 - 100*sum(W_L(index,index,:),[1 2 3])/sum(W_L(index,:,:),[1 2 3]);
    end
end

%% Opening up groups
function posRate = groupOpen(sim_span, X0, psi_t, ZIP, group, STYLE)

    alpha_vec = [logspace(-8,-5,30)];
    K_rep = 10;
    dhetero = [0 5 10 15 20 25];   
    
    for group_index = 1:length(group.name)
        for i = 1:length(alpha_vec)
            alpha = alpha_vec(i);
            [W_L heterophily(i,:)] = LaplaceSmoothing_group(ZIP.W,alpha,group_index,group);
            hetero{group_index}(i) = heterophily(i,group_index);

            
            casePred = [];
            for k = 1:K_rep
                casePred(:,:,k) = SimNet(sim_span, X0, W_L, psi_t, ZIP, group);
            end
            casePred_mean = mean(casePred,3);   % mean exposure across runs
            for g = 1:length(group.name)
                group_pos(g,:)  = sum(casePred_mean(group.ZIPs{g},:),1);
            end
            group_pos(length(group.name)+1,:)  = sum(casePred_mean(:,:),1);
            group_posRate{group_index}(1:length(group.name),i)  = [100*sum(group_pos(1:length(group.name),sim_span),2)./mean(group.pop(1:length(group.name),sim_span),2)]';
            group_posRate{group_index}(length(group.name)+1,i)  = [100*sum(group_pos(length(group.name)+1,sim_span),2)/mean(sum(group.pop(:,sim_span),1),2)]';
        end
    end
    
    for group_index = 1:length(group.name)
        for g = 1:length(group.name)
            for i = 1:length(dhetero)
                posRate{group_index}(g,i) = interp1(hetero{group_index},group_posRate{group_index}(g,:),hetero{group_index}(1)+dhetero(i));
            end
        end
    end
end


%% Isolate groups by removing all incoming and outgoing edges and rescaling the rest
function W_new = IsolateGroups(W,group_ZIPs,groups_to_isolate)   
    W_new = W;
    for t = 1:size(W,3)
        for g = 1:length(groups_to_isolate)
            index = group_ZIPs{groups_to_isolate(g)};
            
            % reassign all outgoing edges
            for i = 1:length(index)               
                dW  = W(index(i),:,t);
                dW2 = W(index(i),index,t);
                fac = sum(dW)/sum(dW2);
                
                W_new(index(i),:,t) = 0;
                W_new(index(i),index,t) = fac*dW2;
            end
            
            % remove all incoming edges not from the group
            for i = 1:length(index)               
                w  = W_new(index,index(i),t);
                W_new(:,index(i),t) = 0;
                W_new(index,index(i),t) = w;
            end

        end
        for j = 1:size(W,2)
            fac = sum(W(j,:,t),'omitnan')/sum(W_new(j,:,t),'omitnan');
            W_new(j,:,t) = fac*W_new(j,:,t);
        end
    end
end


%% Isolate ZIPS by removing all incoming and outgoing edges and rescaling the rest
function W_new = IsolateNodes(W,group_ZIPs,groups_to_isolate)
    W_new = W;
    for t = 1:size(W,3)
        for g = 1:length(groups_to_isolate)
            index = group_ZIPs{groups_to_isolate(g)};
            % reassign all outgoing edges
            for i = 1:length(index) 
                W_new(index(i),:,t) = 0;
                W_new(index(i),index(i),t) = sum(W(index(i),:,t),2);
            end
            
            % remove all incoming edges
            for i = 1:length(index) 
                w = W_new(index(i),index(i),t);
                W_new(:,index(i),t) = 0;
                W_new(index(i),index(i),t) = w;
            end
        end
        for j = 1:size(W,2)
            fac = sum(W(j,:,t),'omitnan')/sum(W_new(j,:,t),'omitnan');
            W_new(j,:,t) = fac*W_new(j,:,t);
        end
        
    end
end

%% Isolate groups and nodes
function [group_rate heterophily] = Isolation(sim_span, X0, W, psi_t, ZIP, group, STYLE)

    K_rep = 10;

    % Without any modifications
    for k = 1:K_rep
        casePred(:,:,k) = SimNet(sim_span, X0, W, psi_t, ZIP, group);
    end
    casePred = mean(casePred,3);
    for g = 1:length(group.name)
       group_rate(1:4,1,g) = 100*sum(casePred(group.ZIPs{g},:),[1 2])/group.pop_mean(g);
    end
    group_rate(1:4,1,length(group.name)+1) = 100*sum(casePred,[1 2])/sum(group.pop_mean);


    for g = 1:length(group.name)
        index = group.ZIPs{g};
        heterophily.original(1,g) = 100 - 100*sum(W(index,index,:),[1 2 3],'omitnan')/sum(W(index,:,:),[1 2 3],'omitnan');
    end
    

    for i = 1:length(group.name)
        % Isolate group i
        W_new = IsolateGroups(ZIP.W,group.ZIPs,[i]);
        for g = 1:length(group.name)
            index = group.ZIPs{g};
            heterophily.groups(i,g) = 100 - 100*sum(W_new(index,index,:),[1 2 3],'omitnan')/sum(W_new(index,:,:),[1 2 3],'omitnan');
        end

        for k = 1:K_rep
            casePred(:,:,k) = SimNet(sim_span, X0, W_new, psi_t, ZIP, group);
        end
        casePred = mean(casePred,3);
        for g = 1:length(group.name)
           group_rate(i,2,g) = 100*sum(casePred(group.ZIPs{g},:),[1 2])/group.pop_mean(g); 
        end
        group_rate(i,2,length(group.name)+1) = 100*sum(casePred,[1 2])/sum(group.pop_mean);

        % Isolate nodes in group i
        W_new = IsolateNodes(ZIP.W,group.ZIPs,[i]);
        for g = 1:length(group.name)
            index = group.ZIPs{g};
            heterophily.nodes(i,g) = 100 - 100*sum(W_new(index,index,:),[1 2 3],'omitnan')/sum(W_new(index,:,:),[1 2 3],'omitnan');
        end
        for k = 1:K_rep
            casePred(:,:,k) = SimNet(sim_span, X0, W_new, psi_t, ZIP, group);
        end
        casePred = mean(casePred,3);
        for g = 1:length(group.name)
           group_rate(i,3,g) = 100*sum(casePred(group.ZIPs{g},:),[1 2])/group.pop_mean(g); 
        end   
        group_rate(i,3,length(group.name)+1) = 100*sum(casePred,[1 2])/sum(group.pop_mean);
    end
end

%% Plot results of changes in mobility patterns (Laplace and isolation)
function PlotStructure(sim_span, X0, psi, ZIP, group, STYLE)
    [parentdir,~,~]=fileparts(pwd);

    posRate = groupOpen(sim_span, X0, psi, ZIP, group, STYLE);
    [group_rate heterophily] = Isolation(sim_span, X0, ZIP.W, psi, ZIP, group, STYLE);
    MarkerSize = 150;
    
    % posRate coming from opening via Laplace smoothing
    % posRate{group_index}(g,i) means the following
    % case rate of group g when smoothing group_index group such that its homophily
    % decreases with (i-1)*5 %. (first entry is baseline)
    
    % heterophily is coming from isolation
    % heterophily.orginial is baseline for each group
    % heterophily.groups(i,j) is for group j when group i is isolated
    % heterophily.nodes(i,j) is for group j when nodes in group are isolated
    
    % group_rate is coming from isolation
    % group_rate(g,i,j) is case rate of group j (5 is city)
    % when group g is manipulated such that 
    % i = 1: baseline, i = 2: groups isolated, i = 3: nodes isolated
    
    figure('Position',[0 0 250 250])
    tlt = tiledlayout(1, 1);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';
    
    
    for g = 1:length(group.name)
        hold on
    
        % baseline
        scatter(100-heterophily.original(1,g),group_rate(g,1,g),MarkerSize,'MarkerFaceColor',STYLE.colors.dark{g},'MarkerEdgeColor','k')
        
        scatter(100-heterophily.groups(g,g),group_rate(g,2,g),MarkerSize,'d','MarkerFaceColor',STYLE.colors.light{g},'MarkerEdgeColor','k')
        scatter(100-heterophily.nodes(g,g),group_rate(g,3,g),MarkerSize,'p','MarkerFaceColor',STYLE.colors.light{g},'MarkerEdgeColor','k')
    
        for i = 2:6
            scatter(100-heterophily.original(1,g)-(i-1)*5,posRate{g}(g,i),MarkerSize,'MarkerFaceColor',STYLE.colors.light{g},'MarkerEdgeColor','none')
        end    
        grid on
        xlim([50 100]), ylim([0 100])
        xlabel('homophily (%)'), ylabel('case rate (%)')
        
    end
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize) 
       
    set(gcf,'renderer','painters');
    print(gcf,[parentdir '/figures/structure.eps'],'-depsc')
end











%% Eliminate shift in homophily
function [W_sameHomo homo_before homo_after sameHomo_before sameHomo_after] = Homophily_matching(ZIP,group,before,after,groupSelect)
    N_nodes = size(ZIP.W,2);

    alpha_vec = logspace(-8,-6,30);

    for g = 1:length(group.ZIPs)

        days = before.days;
        homo_before(g) = mean(squeeze(sum(ZIP.W(group.ZIPs{g},group.ZIPs{g},days),[1 2],'omitnan'))'./group.trip(g,days),'omitnan');

        days = after.days;
        homo_after(g) = mean(squeeze(sum(ZIP.W(group.ZIPs{g},group.ZIPs{g},days),[1 2],'omitnan'))'./group.trip(g,days),'omitnan');

        for i = 1:length(alpha_vec)
            alpha = alpha_vec(i);
            WL = [];
            W0 = ZIP.W;

                
            for t = after.days
                WL(:,:,t) = W0(:,:,t);
                index = group.ZIPs{g};
                w = W0(index,:,t);
                d = repmat(sum(w,2),1,N_nodes);
                w_new = (w+alpha*d.^2/N_nodes)./(1+alpha*d);
                WL(index,:,t) = w_new;
            end

            days = after.days;
            homo_hypo(g,i) = mean(squeeze(sum(WL(group.ZIPs{g},group.ZIPs{g},days),[1 2],'omitnan'))'./group.trip(g,days),'omitnan');
        end

        [~ , alpha_index] = min(abs(homo_hypo(g,:) - homo_before(g))); 
        alpha_opt(g) = alpha_vec(alpha_index);
    end


    W_sameHomo = ZIP.W;
    for g = groupSelect
        alpha = alpha_opt(g);
        WL = [];
        W0 = ZIP.W;
        for t = after.days
            WL(:,:,t) = W0(:,:,t);
            index = group.ZIPs{g};
            w = W0(index,:,t);
            d = repmat(sum(w,2),1,N_nodes);
            w_new = (w+alpha*d.^2/N_nodes)./(1+alpha*d);
            W_sameHomo(index,:,t) = w_new;
        end
    end


    for g = 1:length(group.ZIPs)

        days = before.days;
        sameHomo_before(g) = mean(squeeze(sum(W_sameHomo(group.ZIPs{g},group.ZIPs{g},days),[1 2],'omitnan'))'./group.trip(g,days),'omitnan');

        days = after.days;
        sameHomo_after(g) = mean(squeeze(sum(W_sameHomo(group.ZIPs{g},group.ZIPs{g},days),[1 2],'omitnan'))'./group.trip(g,days),'omitnan');
    end
end




%% Eliminate homophily shift (upon March 1)
function [difference inequality homo_before homo_after] = EliminateHomophilyShift(ZIP,group,sim_span,X0,groupSelect)
    before.days = [1:75];
    after.days  = [76:360]; 

    [W_sameHomo homo_before homo_after sameHomo_before sameHomo_after] = Homophily_matching(ZIP,group,before,after,groupSelect);

    for k = 1:100
        casePred = SimNet(sim_span, X0, ZIP.W, ZIP.psi, ZIP, group);
        for g = 1:4
            group_caseRate(k,g) = 100*sum(casePred(group.ZIPs{g},:),[1 2])/group.pop_mean(g);
        end
        group_caseRate(k,5) = 100*sum(casePred,[1 2])/sum(group.pop_mean);

    end
    groupRate.mean = mean(group_caseRate);
    groupRate.std  = std(group_caseRate);
    
    
    for k = 1:100
        casePred = SimNet(sim_span, X0, W_sameHomo, ZIP.psi, ZIP, group);
        for g = 1:4
            group_caseRate_sameHomo(k,g) = 100*sum(casePred(group.ZIPs{g},:),[1 2])/group.pop_mean(g);
        end
        group_caseRate_sameHomo(k,5) = 100*sum(casePred,[1 2])/sum(group.pop_mean);

    end
    groupRate_sameHomo.mean = mean(group_caseRate_sameHomo);
    groupRate_sameHomo.std  = std(group_caseRate_sameHomo);
    
    difference.mean = mean(group_caseRate_sameHomo - group_caseRate);
    difference.std  = std(group_caseRate_sameHomo - group_caseRate);
        
    inequality.mean = mean((max(group_caseRate_sameHomo,[],2)-min(group_caseRate_sameHomo,[],2))-(max(group_caseRate,[],2)-min(group_caseRate,[],2)));
    inequality.std  = std((max(group_caseRate_sameHomo,[],2)-min(group_caseRate_sameHomo,[],2))-(max(group_caseRate,[],2)-min(group_caseRate,[],2)));

    groupRate.mean
    groupRate_sameHomo.mean
end



%% hex2rgb
function [ rgb ] = hex2rgb(hex,range)

    %% Input checks:
    assert(nargin>0&nargin<3,'hex2rgb function must have one or two inputs.') 
    if nargin==2
        assert(isscalar(range)==1,'Range must be a scalar, either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
    end
    %% Tweak inputs if necessary: 
    if iscell(hex)
        assert(isvector(hex)==1,'Unexpected dimensions of input hex values.')
        
        % In case cell array elements are separated by a comma instead of a
        % semicolon, reshape hex:
        if isrow(hex)
            hex = hex'; 
        end
        
        % If input is cell, convert to matrix: 
        hex = cell2mat(hex);
    end
    if strcmpi(hex(1,1),'#')
        hex(:,1) = [];
    end
    if nargin == 1
        range = 1; 
    end
    %% Convert from hex to rgb: 
    switch range
        case 1
            rgb = reshape(sscanf(hex.','%2x'),3,[]).'/255;
        case {255,256}
            rgb = reshape(sscanf(hex.','%2x'),3,[]).';
        
        otherwise
            error('Range must be either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
    end
end