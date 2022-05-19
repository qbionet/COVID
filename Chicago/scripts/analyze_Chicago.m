% Analyze Chicago
 
clear all
close all
clc

[parentdir,~,~]=fileparts(pwd);

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

STATS = load([parentdir '/statistics/stats_combined.mat']);
STATS_noComb = load([parentdir '/statistics/stats_original.mat']);


ZIP   = STATS.ZIP;
group = STATS.group;

ZIP_noComb   = STATS_noComb.ZIP;
group_noComb = STATS_noComb.group;

[value index] = max( ZIP.caseEst~=0, [], 2);
t_start = min(index - 1);
t_end   = 360;
sim_span = t_start:t_end;

% Initial condition 
p0 = 1e-3;             
p0 = [1-p0 p0 0 0];
X0 = ZIP.pop(:,t_start).*p0;
X0_noComb = ZIP_noComb.pop(:,t_start).*p0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MODEL FIT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Original ZIPs
[psi_noComb error_noComb] = Find_psi(sim_span, X0_noComb, ZIP_noComb, group_noComb, STYLE);
[ZIP_noComb, group_noComb] = SpreadStats(sim_span, ZIP_noComb, group_noComb, psi_noComb);
SelfTripFrac_posBars(ZIP_noComb.W, 'N', ZIP_noComb, group_noComb, STYLE);
PlotFit(sim_span, X0_noComb, ZIP_noComb, group_noComb, 'N', STYLE);


%% Combine ZIPs
CombineZIPs(STATS_noComb, STATS, STYLE)

%% Combined ZIPs
[psi error] = Find_psi(sim_span, X0, ZIP, group, STYLE);
[ZIP, group] = SpreadStats(sim_span, ZIP, group, psi);
SelfTripFrac_posBars(ZIP.W, 'Y', ZIP, group, STYLE);
error = PlotFit(sim_span, X0, ZIP, group, 'Y', STYLE);
corr_coeff = betaError(sim_span, X0, ZIP.W, psi, ZIP, group, STYLE);

[mean(error(group.ZIPs{1})) mean(error(group.ZIPs{2})) mean(error(group.ZIPs{3})) mean(error(group.ZIPs{4}))]
[std(error(group.ZIPs{1})) std(error(group.ZIPs{2})) std(error(group.ZIPs{3})) std(error(group.ZIPs{4}))]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  SOURCE OF EXPOSURE  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Where do exposures happen
ZIP_stats_noComb = ExposureSource(sim_span, X0_noComb, ZIP_noComb, group_noComb, psi_noComb, STYLE);
[ZIP_stats group_stats CS] = ExposureSource(sim_span, X0, ZIP, group, psi, STYLE);
PlotHomophilyExposure(sim_span, ZIP, group, CS, STYLE);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Matching  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot beta trajectories
plotBeta(ZIP,group,STYLE);

%% Plot initial scatters
Scatter_initial_and_groupBars(sim_span,X0,ZIP,group,STYLE);

% Violin plot for group level matching
PlotViolin(sim_span, X0, ZIP, group, STYLE);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Structure  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Effects of structural changes (Laplace and isolation)
PlotStructure(sim_span, X0, psi, ZIP, group, STYLE);

%% What would have happened if homophily had not changed compared to before March 1
shift_mean = [];
shift_std  = [];
base_group_posRate = [group.posRate group.posRate*group.pop_mean/sum(group.pop_mean)];
base_inequality = max(base_group_posRate)-min(base_group_posRate);

groupSelect = [1];
[difference inequality homo_before homo_after] = EliminateHomophilyShift(ZIP,group,sim_span,X0,groupSelect);
shift_mean = [shift_mean ; difference.mean inequality.mean];
shift_std  = [shift_std ; difference.std inequality.std];

groupSelect = [2];
[difference inequality homo_before homo_after] = EliminateHomophilyShift(ZIP,group,sim_span,X0,groupSelect);
shift_mean = [shift_mean ; difference.mean inequality.mean];
shift_std  = [shift_std ; difference.std inequality.std];

groupSelect = [3];
[difference inequality homo_before homo_after] = EliminateHomophilyShift(ZIP,group,sim_span,X0,groupSelect);
shift_mean = [shift_mean ; difference.mean inequality.mean];
shift_std  = [shift_std ; difference.std inequality.std];

groupSelect = [4];
[difference inequality homo_before homo_after] = EliminateHomophilyShift(ZIP,group,sim_span,X0,groupSelect);
shift_mean = [shift_mean ; difference.mean inequality.mean];
shift_std  = [shift_std ; difference.std inequality.std];

groupSelect = [1:4];
[difference inequality homo_before homo_after] = EliminateHomophilyShift(ZIP,group,sim_span,X0,groupSelect);
shift_mean = [shift_mean ; difference.mean inequality.mean];
shift_std  = [shift_std ; difference.std inequality.std];

homo_beforeafter = [100*homo_before ; 100*homo_after]






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Tradeoff  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Change the beta randomly of each group 
corr_rb = groupRateBeta(sim_span, X0, ZIP, group, STYLE);


%% Census scatter
corr_cen = CensusScatter_horizontal(ZIP,group,STYLE);

%% Trade-offs and census changes
CensusTripBars(ZIP,group,corr_rb,corr_cen,STYLE);
PlotTradeOffs(ZIP,group,corr_rb,corr_cen,STYLE);












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MODEL FIT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Calculate psi
function [psi_t error_t] = Find_psi(sim_span, X0,  ZIP, group, STYLE)
    
    delta_E = 4;            % mean latency period (4 days)
    delta_I = 3.5;          % mean infectious period (3.5 days)
    
    t_start = sim_span(1);
    t_end   = sim_span(end); 

    nodes = size(X0,1);
    % initial condition
    X(1:nodes,1:4,t_start) = X0;
    
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

        % calculate psi_t for each node
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
        for g = 1:4
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
        
        for g1 = 1:4
            index1 = group.ZIPs{g1};
            for g2 = 1:4
                index2 = group.ZIPs{g2};
                risk.group(g1,g2,t) = sum(W_t(index1,index2)*lambda_t(index2))/sum(W_t(index1,:)*lambda_t(:));
                
                
            end
        end
        for i = 1:nodes
            for g = 1:4
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
    for g = 1:4
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



%% Self trip vs average other trips and positivity bars
function SelfTripFrac_posBars(W, combined, ZIP, group, STYLE)

    [parentdir,~,~]=fileparts(pwd);

    figure('Position',[0 0 1000 150])
    tlt = tiledlayout(1, 4);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';

    nexttile
    for t = 1:size(W,3)
       dW = triu(W(:,:,t),1)+tril(W(:,:,t),-1);
       row_mean = sum(dW,2)/(size(dW,2)-1);
       w_frac(:,t) = diag(W(:,:,t))./row_mean;
    end

    [~,edges] = histcounts(log10(w_frac(:)),250);
    histogram(w_frac(:),10.^edges,'EdgeColor','none','FaceColor',hex2rgb(STYLE.colors.dark{5}))
    set(gca,'XScale','log')
    grid on
    xlabel('within-ZIP trips per average between-ZIP trip')
    ylabel('#')
    xlim([1e0 1e4])
    xticks([1e0 1e1 1e2 1e3 1e4])
    
    nexttile
    [value, index] = sort(ZIP.posRate,'descend');

    b = bar(value,'FaceColor','flat','EdgeColor','none');
    for i = 1:length(index)
        for g = 1:4
            member(g) = ismember(index(i),group.ZIPs{g});
        end
        b.CData(i,:) = hex2rgb(STYLE.colors.dark{find(member==1)}); 
    end
    grid on
    xlim([0 length(index)+1])
    ylim([0 100])
    xticks([])
    xlabel('ZIPs')
    ylabel('case rate (%)')
    
    
    
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize) 
    set(gcf, 'Renderer', 'painters');
    switch combined
        case 'N'
            print(gcf,[parentdir '/figures/fit/selfTripFrac_noComb.eps'],'-depsc')
        case 'Y'
            print(gcf,[parentdir '/figures/fit/selfTripFrac.eps'],'-depsc')
    end
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
    tlt = tiledlayout(2, 4);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';   
    for g = 1:4
        nexttile
        plot(sim_span, squeeze(sum(casePred_sims(group.ZIPs{g},sim_span,:),1)),'o','MarkerFaceColor',STYLE.colors.light{g},'MarkerEdgeColor','none','MarkerSize',4)
        hold on
        plot(sim_span, sum(casePred(group.ZIPs{g},sim_span),1),'Color',STYLE.colors.dark{g},'LineWidth',STYLE.thick)
        plot(sim_span, group.caseEst(g,sim_span),'k','LineWidth',STYLE.thin)
        xticks(STYLE.XTick_months)
        xticklabels(STYLE.XTickLabel_months)
        xlim(STYLE.XLim_year)
        ylim([0 6000])
        ylabel('case rate (#)')
        title(group.name{g})
        xtickangle(0)
        grid on
    end
    
    
    
    nexttile([1 2])
    hold on
    for g = 1:4
        index = group.ZIPs{g};
        s = scatter(1e-4*meanPop(group.ZIPs{g}),error(group.ZIPs{g}),50,'filled');

        s.MarkerEdgeColor = 'none';
        s.MarkerFaceColor = STYLE.colors.dark{g};
    end
    grid on
    xlim([0 12])
    ylim([-10 60])
    xlabel('population (10k)')
    ylabel('error (%)')
    
    
    nexttile([1 2])
    hold on
    for g = 1:4
        index = group.ZIPs{g};
        s = scatter(100-ZIP.hetero(group.ZIPs{g}),error(group.ZIPs{g}),50,'filled');

        s.MarkerEdgeColor = 'none';
        s.MarkerFaceColor = STYLE.colors.dark{g};
    end
    grid on
    xlim([0 100])
    ylim([-10 60])
    xlabel('homophily (%)')
    ylabel('error (%)')
    
    set(gcf,'renderer','painters');
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize) 
    switch combined
        case 'N'
            print(gcf,[parentdir '/figures/fit/fit_noComb.eps'],'-depsc')
        case 'Y'
            print(gcf,[parentdir '/figures/fit/fit.eps'],'-depsc')
    end


    figure('Position',[0 0 600 250])
    tlt = tiledlayout(2, 2);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';   
    for g = 1:4
        nexttile
        plot(sim_span, squeeze(sum(casePred_sims(group.ZIPs{g},sim_span,:),1)),'o','MarkerFaceColor',STYLE.colors.light{g},'MarkerEdgeColor','none','MarkerSize',4)
        hold on
        plot(sim_span, sum(casePred(group.ZIPs{g},sim_span),1),'Color',STYLE.colors.dark{g},'LineWidth',STYLE.thick)
        plot(sim_span, group.caseEst(g,sim_span),'k','LineWidth',STYLE.thin)
        xticks(STYLE.XTick_months)
        xticklabels(STYLE.XTickLabel_months)
        xlim(STYLE.XLim_year)
        ylim([0 6000])
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
            print(gcf,[parentdir '/figures/fit/fit_compact_noComb.eps'],'-depsc')
        case 'Y'
            print(gcf,[parentdir '/figures/fit/fit_compact.eps'],'-depsc')
    end
end







%% Combine ZIPs
function CombineZIPs(STATS_noComb, STATS, STYLE)
    [parentdir,~,~]=fileparts(pwd);

    figure('Position',[0 0 1000 150])
    tlt = tiledlayout(1, 4);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';



    nexttile
    [area_noComb index_noComb] = sort(STATS_noComb.ZIP.area(STATS_noComb.group.ZIPs{3}),'descend');
    [area_comb index_comb]     = sort(STATS.ZIP.area(STATS.group.ZIPs{3}),'descend');
    area_comb                  = [area_comb ; zeros(length(area_noComb)-length(area_comb),1)];
    b = bar([area_noComb area_comb]','FaceColor','flat','EdgeColor','none');
    for i = 1:length(b)
        b(i).CData = repmat(hex2rgb(STYLE.colors.dark{5}),2,1); 
    end
    name_noComb = STATS_noComb.ZIP.name(STATS_noComb.group.ZIPs{3});
    name_noComb = name_noComb(index_noComb);
    for i = 1:length(STATS.ZIPs_to_combine)
        index = find(name_noComb == STATS.ZIPs_to_combine(i));
        b(index).CData(1,:) = hex2rgb(STYLE.colors.dark{3});
    end
    name_comb = STATS.ZIP.name(STATS.group.ZIPs{3});
    name_comb = name_comb(index_comb);
    index = find(name_comb == STATS.ZIPs_to_combine(1));
    b(index).CData(2,:) = hex2rgb(STYLE.colors.dark{3});
    xticklabels({'original' 'combined'})
    ylabel('area (sq mile)')
    ylim([0 10])
    grid on

    nexttile
    [pop_noComb index_noComb] = sort(STATS_noComb.ZIP.pop(STATS_noComb.group.ZIPs{3}),'descend');
    [pop_comb index_comb]     = sort(STATS.ZIP.pop(STATS.group.ZIPs{3}),'descend');
    pop_comb                  = [pop_comb ; zeros(length(pop_noComb)-length(pop_comb),1)];
    b = bar(1e-4*[pop_noComb pop_comb]','FaceColor','flat','EdgeColor','none');
    for i = 1:length(b)
        b(i).CData = repmat(hex2rgb(STYLE.colors.dark{5}),2,1); 
    end
    name_noComb = STATS_noComb.ZIP.name(STATS_noComb.group.ZIPs{3});
    name_noComb = name_noComb(index_noComb);
    for i = 1:length(STATS.ZIPs_to_combine)
        index = find(name_noComb == STATS.ZIPs_to_combine(i));
        b(index).CData(1,:) = hex2rgb(STYLE.colors.dark{3});
    end
    name_comb = STATS.ZIP.name(STATS.group.ZIPs{3});
    name_comb = name_comb(index_comb);
    index = find(name_comb == STATS.ZIPs_to_combine(1));
    b(index).CData(2,:) = hex2rgb(STYLE.colors.dark{3});
    xticklabels({'original' 'combined'})
    ylabel('population (10k)')
    grid on


    nexttile
    hold on
    plot(STATS_noComb.ZIP.tripRate(STATS.comb,:)','LineWidth',STYLE.thin,'Color',STYLE.colors.light{3})
    plot(STATS.ZIP.tripRate(STATS.comb(1),:)','LineWidth',STYLE.thick,'Color',STYLE.colors.dark{3})
    plot(STATS.group.tripRate(3,:),'k')
    grid on
    xticks(STYLE.XTick_months)
    xticklabels(STYLE.XTickLabel_months)
    xlim(STYLE.XLim_year)
    ylabel('trip rate')
    ylim([0 20])
    xtickangle(0)

    nexttile
    hold on
    plot(100*[STATS_noComb.ZIP.caseEst(STATS.comb,:)./STATS_noComb.ZIP.pop(STATS.comb,:)]','LineWidth',STYLE.thin,'Color',STYLE.colors.light{3})
    plot(100*STATS.ZIP.caseEst(STATS.comb(1),:)./STATS.ZIP.pop(STATS.comb(1),:),'LineWidth',STYLE.thick,'Color',STYLE.colors.dark{3})
    plot(100*STATS.group.caseEst(3,:)./STATS.group.pop(3,:),'k')
    grid on
    xticks(STYLE.XTick_months)
    xticklabels(STYLE.XTickLabel_months)
    xlim(STYLE.XLim_year)
    ylabel('case rate (%)')
    ylim([0 1])
    xtickangle(0)

    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize)
    set(gcf, 'Renderer', 'painters');
    print(gcf,[parentdir '/figures/fit/combineZIPs_stats.eps'],'-depsc')
end


%% beta error 
function corr_coeff = betaError(sim_span, X0, W, psi_t, ZIP, group, STYLE)
    [parentdir,~,~]=fileparts(pwd);

    K_rep = 10;

    for k = 1:K_rep
        [casePred fractions caseSource risk X] = SimNet(sim_span, X0, W, psi_t, ZIP, group);
        [row,column,depth] = size(X);
        X(:,:,1) = [];
        
        for i = 1:size(W,1)
            S = squeeze(X(i,1,:));
            I = squeeze(X(i,3,:));
            N = squeeze(ZIP.pop(i,:))';
            C = squeeze(ZIP.caseEst(i,:))';
            betaSim_rep(i,:,k) = C.*N./(S.*I);
            
        end
        
    end
    betaSim = squeeze(mean(betaSim_rep,3));
    betaSim(isnan(betaSim)) = 0;
    betaSim(isinf(betaSim)) = 0;
        
    
    for i = 1:size(W,1)
        beta(i,:) = psi_t(i,:).*ZIP.popDens(i,:).*ZIP.selfRate(i,:).^2;
    
    end

    beta_mean    = mean(beta(:,sim_span),2,'omitnan');
    betaSim_mean = mean(betaSim(:,sim_span),2,'omitnan');
    betaDiff = 100*abs(betaSim-beta)./betaSim;    
    
    diffMean = mean(betaDiff(:,sim_span),2,'omitnan');
        
    figure('Position',[0 0 250 250])
    tlt = tiledlayout(1, 1);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';
    
    hold on
    plot(beta_mean,betaSim_mean,'bo','MarkerFaceColor','b')
    plot([0 1],[0 1],'k:','LineWidth',3)
    xlim([0 1])
    ylim([0 1])
    axis square
    grid on
    xticks([0:0.2:1])
    yticks([0:0.2:1])
    
    
    xlabel('$\beta$','Interpreter','latex')
    ylabel('$\bar \beta$','Interpreter','latex')
    
    delta = linspace(0,1,1e3);
    for i = 1:length(delta)
       error(i) = sum((beta_mean*(1+delta(i)) - betaSim_mean).^2); 
    end
    
    C = corrcoef([beta_mean betaSim_mean]);
    corr_coeff = C(1,2);
 
    % best fit
    [value, index] = min(error);
    [delta(index) sqrt(1+delta(index))];

    betai = linspace(0,1,1e2);
    plot(betai,betai*(1+delta(index)),'m','LineWidth',3)
    
    set(gcf,'renderer','painters');
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize) 
    print(gcf,[parentdir '/figures/fit/betaDiff.eps'],'-depsc')
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  SOURCE OF EXPOSURE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Trip risk
function [ZIP_stats group_stats CS] = ExposureSource(sim_span, X0, ZIP, group, psi_t, STYLE)
    [parentdir,~,~]=fileparts(pwd);

    K_rep = 100;

    figure('Position',[0 0 1000 250])
    tlt = tiledlayout(1, 4);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';  


    % Compute where exposures happen
    for k = 1:K_rep
        [casePred fractions(k) cs risk] = SimNet(sim_span, X0, ZIP.W, psi_t, ZIP, group);


        % cs is casesource: n x n x t, each row contains the # of exposures
        % from all nodes at time t

        CS{k} = cs;
        selfRisk(:,:,k) = risk.self;
        for g = 1:4
            groupRisk(g,:,k) = squeeze(risk.group(g,g,:));
        end
        
        for g = 1:4
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
    for g = 1:4
        index = group.ZIPs{g};       
        w_self  = squeeze(sum(ZIP.W(index,index,:),[1 2]));
        w_total = squeeze(sum(ZIP.W(index,:,:),[1 2]));
        groupTrips(g,:) = w_self./w_total;
    end
    
    % Self trips
    nexttile
    hold on
    for g = 1:4
        index = group.ZIPs{g};
        selfTrips_groupMean(g,:) = mean(selfTrips(index,:),1);
    end
    for g = 1:4
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
    for g = 1:4
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
    for g = 1:4
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
    for g = 1:4
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
    [row,column] = size(stats);
    for i = 1:column
        stats(row+1,i) = group.pop_mean'*stats(1:row,i)/sum(group.pop_mean);
    end
    
    for g = 1:4
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
    print(gcf,[parentdir '/figures/source/exposureSource.eps'],'-depsc')
    
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
    print(gcf,[parentdir '/figures/source/homophily_exposure_circles.eps'],'-depsc','-r1200')

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Matching  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot beta for ZIPs and groups
function plotBeta(ZIP,group,STYLE)
    [parentdir,~,~]=fileparts(pwd);

    figure('Position',[0 0 1000 250])
    tlt = tiledlayout(1, 4);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';
    
    beta_mean = mean(ZIP.beta,1);
    beta_std  = std(ZIP.beta,1);
    
    for i = 1:size(ZIP.W,1)
        signal(i,:) = ZIP.beta(i,:)./beta_mean;
    end
    for t = 1:size(ZIP.W,3)
        [value_ascend, index_ascend] = sort(signal(:,t),'ascend');
        [value_descend, index_descend] = sort(signal(:,t),'descend');
        LB(t) = value_ascend(5);
        UB(t) = value_descend(5);
    end
    
    nexttile
    hold on
    plot(LB,'LineWidth',STYLE.medium,'Color',STYLE.colors.dark{5})
    plot(UB,'LineWidth',STYLE.medium,'Color',STYLE.colors.dark{5})
    axis square
    grid on    
    ylim([0 2])
    xlim(STYLE.XLim_year)
    xticks(STYLE.XTick_months)
    xticklabels(STYLE.XTickLabel_months)
    xtickangle(0)
    ylabel('\beta_i / \beta_{mean}')
    
    set(gcf,'renderer','painters');
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize) 
    print(gcf,[parentdir '/figures/matching/betaRange.eps'],'-depsc')
end


%% Plot initial scatters and group bars
function [R P] = Scatter_initial_and_groupBars(sim_span,X0,ZIP,group,STYLE)
    [parentdir,~,~]=fileparts(pwd);

    figure('Position',[0 0 1000 350])
    tlt = tiledlayout(2, 4);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';

    nexttile
    x_ZIP   = ZIP.psi;
    y_ZIP   = ZIP.posRate;
    s_ZIP   = ZIP.pop;
    x_group = group.psi;
    y_group = group.posRate;
    s_group = group.pop;
    PlotScatter(sim_span,x_ZIP,y_ZIP,s_ZIP,x_group,y_group,s_group,group.ZIPs,[],STYLE)
    ylim([0 100])
    xlabel('\psi')
    ylabel('case rate (%)')
    [r p] = corrcoef(mean(x_ZIP(:,sim_span),2)', y_ZIP');
    R(1) = r(1,2);
    P(1) = p(1,2);
    
    
    nexttile
    x_ZIP   = ZIP.popDens/1e4;
    y_ZIP   = ZIP.posRate;
    s_ZIP   = ZIP.pop;
    x_group = group.popDens/1e4;
    y_group = group.posRate;
    s_group = group.pop;
    PlotScatter(sim_span,x_ZIP,y_ZIP,s_ZIP,x_group,y_group,s_group,group.ZIPs,[],STYLE)
    ylim([0 100])
    xlabel('population density (10k/mile^2)')
    ylabel('case rate (%)')
    [r p] = corrcoef(mean(x_ZIP(:,sim_span),2)', y_ZIP');
    R(2) = r(1,2);
    P(2) = p(1,2);
    
    nexttile
    x_ZIP   = ZIP.tripRate;
    y_ZIP   = ZIP.posRate;
    s_ZIP   = ZIP.pop;
    x_group = group.tripRate;
    y_group = group.posRate;
    s_group = group.pop;
    PlotScatter(sim_span,x_ZIP,y_ZIP,s_ZIP,x_group,y_group,s_group,group.ZIPs,[],STYLE)
    ylim([0 100])
    xlabel('tripRate')
    ylabel('case rate (%)')
    [r p] = corrcoef(mean(x_ZIP(:,sim_span),2)', y_ZIP');
    R(3) = r(1,2);
    P(3) = p(1,2);
    
    nexttile
    x_ZIP   = ZIP.beta;
    y_ZIP   = ZIP.posRate;
    s_ZIP   = ZIP.pop;
    x_group = group.beta;
    y_group = group.posRate;
    s_group = group.pop;
    PlotScatter(sim_span,x_ZIP,y_ZIP,s_ZIP,x_group,y_group,s_group,group.ZIPs,[],STYLE)
    ylim([0 100])
    xlabel('\beta')
    ylabel('case rate (%)')
    [r p] = corrcoef(mean(x_ZIP(:,sim_span),2)', y_ZIP');
    R(4) = r(1,2);
    P(4) = p(1,2);


    nexttile
    b = bar(group.posRate/max(group.posRate),'FaceColor','flat','EdgeColor','none');
    for g = 1:4
        b.CData(g,:) = hex2rgb(STYLE.colors.dark{g});
    end
    grid on
    xlim([0.5 4.5])
    ylim([0 1.1])
    xticklabels({'Black' 'Latinx' 'White' 'Mixed'})
    title('case rate')

    nexttile
    b = bar(group.beta_mean/max(group.beta_mean),'FaceColor','flat','EdgeColor','none');
    for g = 1:4
        b.CData(g,:) = hex2rgb(STYLE.colors.dark{g});
    end
    grid on
    xlim([0.5 4.5])
    ylim([0 1.1])
    xticklabels({'Black' 'Latinx' 'White' 'Mixed'})
    title('\beta')

    nexttile
    b = bar(group.beta1_mean/max(group.beta1_mean),'FaceColor','flat','EdgeColor','none');
    for g = 1:4
        b.CData(g,:) = hex2rgb(STYLE.colors.dark{g});
    end
    grid on
    xlim([0.5 4.5])
    ylim([0 1.1])
    xticklabels({'Black' 'Latinx' 'White' 'Mixed'})
    title('vulnerability')

    nexttile
    b = bar(group.beta2_mean/max(group.beta2_mean),'FaceColor','flat','EdgeColor','none');
    for g = 1:4
        b.CData(g,:) = hex2rgb(STYLE.colors.dark{g});
    end
    grid on
    xlim([0.5 4.5])
    ylim([0 1.1])
    xticklabels({'Black' 'Latinx' 'White' 'Mixed'})
    title('(trip rate)^2')

    set(gcf,'renderer','painters');
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize) 
    print(gcf,[parentdir '/figures/matching/initialScatter_groupBars.eps'],'-depsc')
end


%% Plot scatters
function PlotScatter(sim_span,x_ZIP,y_ZIP,s_ZIP,x_group,y_group,s_group,group_ZIPs,match,STYLE)    
    hold on
    for g = 1:4
       for i = 1:length(group_ZIPs{g})
           index = group_ZIPs{g}(i);
           scatter(mean(x_ZIP(index,sim_span)),y_ZIP(index),0.2*log10(mean(s_ZIP(index,sim_span)))^4,hex2rgb(STYLE.colors.light{g}),'filled')
       end
    end
    for g = 1:4
        s = scatter(mean(x_group(g,sim_span)),y_group(g),0.2*log10(mean(s_group(g,sim_span)))^4,'filled');
        s.MarkerEdgeColor = 'k';
        s.MarkerFaceColor = hex2rgb(STYLE.colors.dark{g});
        s.LineWidth = 1;
    end
    grid on
    
    if ~isempty(match)
        xdata = [mean(x_group(match.source,sim_span)) mean(x_group(match.target,sim_span))];
        ydata = [y_group(match.source) y_group(match.target)];
        plot(xdata,ydata,'k','LineWidth',2); 
    end
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
    hypoStats.ZIP.posRate = 100*sum(casePred_mean,2)./mean(ZIP.pop,2);
    for g = 1:4
        group_pos(g,:)  = sum(casePred_mean(group.ZIPs{g},:),1);
    end
    hypoStats.group.posRate  = 100*sum(group_pos,2)./mean(group.pop,2);
    hypoStats.group.posCum   = 100*cumsum(group_pos,2)./mean(group.pop,2);

    % calculate trip rate for ZIPs and groups
    ZIP_trip_hypo = squeeze(sum(W,2));
    hypoStats.ZIP.tripRate = ZIP_trip_hypo./ZIP.pop;
    dW = squeeze(sum(W,2));
    for g = 1:4
        group_trip_hypo(g,:) = sum(dW(group.ZIPs{g},:),1);
    end
    hypoStats.group.tripRate = group_trip_hypo./group.pop;

    
    % calculate pop dens for ZIPS and groups
    for i = 1:length(ZIP.area)
        hypoStats.ZIP.popDens(i,:) = ZIP.pop(i,:)./ZIP.area(i);
    end
    for g = 1:4
        index = group.ZIPs{g};
        hypoStats.group.popDens(g,:)  = group.pop(g,:)./sum(ZIP.area(index));
    end
    
    
    % calculate psi for ZIPS and groups
    hypoStats.ZIP.psi = psi;
    for g = 1:4
        index = group.ZIPs{g};
        hypoStats.group.psi(g,:) = sum(psi(index,:).*ZIP.pop(index,:),1)./sum(ZIP.pop(index,:),1);
    end
    
    % calculate vulnerability for ZIPS and groups
    hypoStats.ZIP.vuln = psi.*hypoStats.ZIP.popDens;
    for g = 1:4
        index = group.ZIPs{g};
        hypoStats.group.vuln(g,:) = sum(hypoStats.ZIP.vuln(index,:).*ZIP.pop(index,:),1)./sum(ZIP.pop(index,:),1);
    end
    
    
    % calculate beta for ZIPs and groups
    for i = 1:length(ZIP.area)
        ZIP_selfRate(i,:) = squeeze(W(i,i,:))'./ZIP.pop(i,:);
        hypoStats.ZIP.beta(i,:) = hypoStats.ZIP.psi(i,:).*hypoStats.ZIP.popDens(i,:).*ZIP_selfRate(i,:).^2;
    end
    for g = 1:4
        index = group.ZIPs{g};
        group_selfRate(g,:) = squeeze(sum(W(index,index,:),[1,2]))'./group.pop(g,:);
        hypoStats.group.beta(g,:) = sum(hypoStats.ZIP.beta(index,:).*ZIP.pop(index,:),1)./sum(ZIP.pop(index,:),1);
    end
    
end



%% Matching scatter plots
function rateDiff = MatchingScatters(sim_span,X0,ZIP,group,match,STYLE)
    [parentdir,~,~]=fileparts(pwd);

    ZIPs_toRescale = group.ZIPs{match.source};
    
    figure('Position',[0 0 1000 500])
    tlt = tiledlayout(3, 4);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';
    
    rateDiff(1) = group.posRate(match.target) - group.posRate(match.source);
      

    %% case count as a function of trip rate
    nexttile
    PlotScatter(sim_span,ZIP.tripRate,ZIP.posRate,ZIP.pop,group.tripRate,group.posRate,group.pop,group.ZIPs,match,STYLE);
    xlim([0 20])
    ylim([0 100])
    xlabel('trip rate')
    ylabel('case count (%)')
    title('baseline')
    
    %% rescale triprate temporal
    W_hypo    = ZIP.W;
    psi_hypo  = ZIP.psi;
    
    for i = 1:length(ZIPs_toRescale)
        index = ZIPs_toRescale(i);
        factor_t = group.tripRate(match.target,:)./ZIP.tripRate(index,:);
        tripRate_old = squeeze(sum(W_hypo(index,:,:),2))./ZIP.pop(index,:);
        for t = 1:length(factor_t)
            W_hypo(index,:,t) = factor_t(t)*W_hypo(index,:,t);
        end
        tripRate_new = squeeze(sum(W_hypo(index,:,:),2))./ZIP.pop(index,:);
        factor_mean = mean(tripRate_old(sim_span))/mean(tripRate_new(sim_span));
        W_hypo(index,:,:) = factor_mean*W_hypo(index,:,:);
    end
    hypoStats = HypoStats(sim_span, X0, W_hypo, psi_hypo, ZIP, group);

    nexttile
    x_ZIP   = hypoStats.ZIP.tripRate;
    y_ZIP   = hypoStats.ZIP.posRate;
    s_ZIP   = ZIP.pop;
    x_group = hypoStats.group.tripRate;
    y_group = hypoStats.group.posRate;
    s_group = group.pop;
    PlotScatter(sim_span,x_ZIP,y_ZIP,s_ZIP,x_group,y_group,s_group,group.ZIPs,match,STYLE)
    xlim([0 20])
    ylim([0 100])
    xlabel('trip rate')
    ylabel('case count (%)')
    title('trajectory matching')
    
    %% rescale triprate mean
    W_hypo    = ZIP.W;
    psi_hypo  = ZIP.psi;
    
    factor_W = group.tripRate_mean(match.target)/group.tripRate_mean(match.source);
    W_hypo(ZIPs_toRescale,:,:) = factor_W*W_hypo(ZIPs_toRescale,:,:);
    
    hypoStats = HypoStats(sim_span, X0, W_hypo, psi_hypo, ZIP, group);

    nexttile
    x_ZIP   = hypoStats.ZIP.tripRate;
    y_ZIP   = hypoStats.ZIP.posRate;
    s_ZIP   = ZIP.pop;
    x_group = hypoStats.group.tripRate;
    y_group = hypoStats.group.posRate;
    s_group = group.pop;
    PlotScatter(sim_span,x_ZIP,y_ZIP,s_ZIP,x_group,y_group,s_group,group.ZIPs,match,STYLE)
    xlim([0 20])
    ylim([0 100])
    xlabel('trip rate')
    ylabel('case count (%)')
    title('mean matching')
      
    
    %% rescale triprate mean and temporal
    W_hypo    = ZIP.W;
    psi_hypo  = ZIP.psi;
    
    for i = 1:length(ZIPs_toRescale)
        index = ZIPs_toRescale(i);
        factor_t = group.tripRate(match.target,:)./ZIP.tripRate(index,:);
        tripRate_old = squeeze(sum(W_hypo(index,:,:),2))./ZIP.pop(index,:);
        for t = 1:length(factor_t)
            W_hypo(index,:,t) = factor_t(t)*W_hypo(index,:,t);
        end
        tripRate_new = squeeze(sum(W_hypo(index,:,:),2))./ZIP.pop(index,:);
        factor_mean = mean(tripRate_old(sim_span))/mean(tripRate_new(sim_span));
        W_hypo(index,:,:) = factor_mean*W_hypo(index,:,:);
    end
    factor_W = group.tripRate_mean(match.target)/group.tripRate_mean(match.source);
    W_hypo(ZIPs_toRescale,:,:) = factor_W*W_hypo(ZIPs_toRescale,:,:);

    hypoStats = HypoStats(sim_span, X0, W_hypo, psi_hypo, ZIP, group);

    nexttile
    x_ZIP   = hypoStats.ZIP.tripRate;
    y_ZIP   = hypoStats.ZIP.posRate;
    s_ZIP   = ZIP.pop;
    x_group = hypoStats.group.tripRate;
    y_group = hypoStats.group.posRate;
    s_group = group.pop;
    PlotScatter(sim_span,x_ZIP,y_ZIP,s_ZIP,x_group,y_group,s_group,group.ZIPs,match,STYLE)
    xlim([0 20])
    ylim([0 100])
    xlabel('trip rate')
    ylabel('case count (%)')
    title('mean and trajectory matching')
    
    rateDiff(2) = hypoStats.group.posRate(match.target) - hypoStats.group.posRate(match.source);

    
    %% case count as a function of vulnerability
    nexttile
    PlotScatter(sim_span,ZIP.vuln,ZIP.posRate,ZIP.pop,group.vuln,group.posRate,group.pop,group.ZIPs,match,STYLE);
    % xlim([1e-8 1e-4])
    ylim([0 100])
    xlabel('vulnerability')
    ylabel('case count (%)')
    title('baseline')
    set(gca,'XScale','log')
    
    
    %% rescale psi temporal
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
        factor_mean = mean(vuln_old(sim_span))/mean(vuln_new(sim_span),'omitnan');
        psi_hypo(index,:) = factor_mean*psi_hypo(index,:);
        psi_hypo(index,isnan(psi_hypo(index,:))) = 0;
    end;
    psi_hypo(isnan(psi_hypo)) = 0;

    hypoStats = HypoStats(sim_span, X0, W_hypo, psi_hypo, ZIP, group);

    nexttile
    x_ZIP   = hypoStats.ZIP.vuln;
    y_ZIP   = hypoStats.ZIP.posRate;
    s_ZIP   = ZIP.pop;
    x_group = hypoStats.group.vuln;
    y_group = hypoStats.group.posRate;
    s_group = group.pop;
    PlotScatter(sim_span,x_ZIP,y_ZIP,s_ZIP,x_group,y_group,s_group,group.ZIPs,match,STYLE)
    ylim([0 100])
    xlabel('vulnerability')
    ylabel('case count (%)')
    title('trajectory matching')
    set(gca,'XScale','log')
    
    %% rescale psi mean    
    W_hypo    = ZIP.W;
    psi_hypo  = ZIP.psi;
    
    factor_psi  = group.vuln_mean(match.target)/group.vuln_mean(match.source);
    psi_hypo(ZIPs_toRescale,:) = factor_psi*psi_hypo(ZIPs_toRescale,:);

    hypoStats = HypoStats(sim_span, X0, W_hypo, psi_hypo, ZIP, group);

    nexttile
    x_ZIP   = hypoStats.ZIP.vuln;
    y_ZIP   = hypoStats.ZIP.posRate;
    s_ZIP   = ZIP.pop;
    x_group = hypoStats.group.vuln;
    y_group = hypoStats.group.posRate;
    s_group = group.pop;
    PlotScatter(sim_span,x_ZIP,y_ZIP,s_ZIP,x_group,y_group,s_group,group.ZIPs,match,STYLE)
    ylim([0 100])
    xlabel('vulnerability')
    ylabel('case count (%)')
    title('mean matching')
    set(gca,'XScale','log')
        
    
    
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
        factor_mean = mean(vuln_old(sim_span))/mean(vuln_new(sim_span),'omitnan');
        psi_hypo(index,:) = factor_mean*psi_hypo(index,:);
    end;
    psi_hypo(isnan(psi_hypo)) = 0;
        
    factor_psi  = group.vuln_mean(match.target)/group.vuln_mean(match.source);
    psi_hypo(ZIPs_toRescale,:) = factor_psi*psi_hypo(ZIPs_toRescale,:);
    
    hypoStats = HypoStats(sim_span, X0, W_hypo, psi_hypo, ZIP, group);

    nexttile
    x_ZIP   = hypoStats.ZIP.vuln;
    y_ZIP   = hypoStats.ZIP.posRate;
    s_ZIP   = ZIP.pop;
    x_group = hypoStats.group.vuln;
    y_group = hypoStats.group.posRate;
    s_group = group.pop;
    PlotScatter(sim_span,x_ZIP,y_ZIP,s_ZIP,x_group,y_group,s_group,group.ZIPs,match,STYLE)
    ylim([0 100])
    xlabel('vulnerability')
    ylabel('case count (%)')
    title('mean and trajectory matching')
    set(gca,'XScale','log')
    
    rateDiff(3) = hypoStats.group.posRate(match.target) - hypoStats.group.posRate(match.source);
       
    
    %% case count as a function of beta
    nexttile
    PlotScatter(sim_span,ZIP.beta,ZIP.posRate,ZIP.pop,group.beta,group.posRate,group.pop,group.ZIPs,match,STYLE);
    xlim([0 0.6])
    xticks([0:0.2:0.6])
    ylim([0 100])
    xlabel('\beta')
    ylabel('case count (%)')
    title('baseline')

    
    %% match beta
    for mt = [0 1 2]
        matchType = mt;
        ZIPs_toRescale = group.ZIPs{match.source};

        W_hypo    = ZIP.W;
        psi_hypo  = ZIP.psi;

        % rescale tripRate temporal
        for i = 1:length(ZIPs_toRescale)
            index = ZIPs_toRescale(i);
            factor_t = group.tripRate(match.target,:)./ZIP.tripRate(index,:);
            tripRate_old = squeeze(sum(W_hypo(index,:,:),2))./ZIP.pop(index,:);
            for t = 1:length(factor_t)
                W_hypo(index,:,t) = factor_t(t)*W_hypo(index,:,t);
            end
            tripRate_new = squeeze(sum(W_hypo(index,:,:),2))./ZIP.pop(index,:);
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
            factor_mean = mean(vuln_old(sim_span))/mean(vuln_new(sim_span),'omitnan');
            psi_hypo(index,:) = factor_mean*psi_hypo(index,:);
            psi_hypo(index,isnan(psi_hypo(index,:))) = 0;
        end;
        psi_hypo(isnan(psi_hypo)) = 0;

        % mean of target beta
        beta_target = group.beta_mean(match.target);

        % modified beta (slight change due to trajectory shifts)
        for i = 1:length(ZIPs_toRescale)
            index = ZIPs_toRescale(i);

            ZIP_popDens(index,:)  = ZIP.pop(index,:)/ZIP.area(index);
            ZIP_selfRate(index,:) = squeeze(W_hypo(index,index,:))'./ZIP.pop(index,:);

            ZIP_beta(index,:) = psi_hypo(index,:).*ZIP_popDens(index,:).*ZIP_selfRate(index,:).^2;
        end
        group_beta(:) = sum(ZIP_beta(ZIPs_toRescale,:).*ZIP.pop(ZIPs_toRescale,:),1)./sum(ZIP.pop(ZIPs_toRescale,:),1);
        beta_source = mean(group_beta(sim_span),2);

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

        
        nexttile
        x_ZIP   = hypoStats.ZIP.beta;
        y_ZIP   = hypoStats.ZIP.posRate;
        s_ZIP   = ZIP.pop;
        x_group = hypoStats.group.beta;
        y_group = hypoStats.group.posRate;
        s_group = group.pop;
        PlotScatter(sim_span,x_ZIP,y_ZIP,s_ZIP,x_group,y_group,s_group,group.ZIPs,match,STYLE)
        xlim([0 0.6])
        xticks([0:0.2:0.6])
        ylim([0 100])
        xlabel('\beta')
        ylabel('case count (%)')
        title(TITLE)
        
        if mt == 2
            rateDiff(4) = hypoStats.group.posRate(match.target) - hypoStats.group.posRate(match.source);
        end
    end
    
    sgtitle(['Matching ' group.name{match.source} ' characteristics to ' group.name{match.target}],'FontWeight','bold')
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize) 
    set(gcf, 'Renderer', 'painters');
    print(gcf,[parentdir '/figures/matching/scatter_' num2str(match.source) num2str(match.target) '.eps'],'-depsc')
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
        H{k} = squeeze(abs(RATEDIFF(:,:,k)));
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
    print(gcf,[parentdir '/figures/matching/matching_violin.eps'],'-depsc','-r1200')
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Structure  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    for g = 1:4
        index = group.ZIPs{g};
        heterophily(g) = 100 - 100*sum(W_L(index,index,:),[1 2 3])/sum(W_L(index,:,:),[1 2 3]);
    end
end

%% Opening up groups
function posRate = groupOpen(sim_span, X0, psi_t, ZIP, group, STYLE)
    [parentdir,~,~]=fileparts(pwd);

    alpha_vec = [logspace(-8,-5,30)];
    K_rep = 10;
    dhetero = [0 5 10 15 20 25];   
    
    for group_index = 1:4
        for i = 1:length(alpha_vec)
            alpha = alpha_vec(i);
            [W_L heterophily(i,:)] = LaplaceSmoothing_group(ZIP.W,alpha,group_index,group);
            hetero{group_index}(i) = heterophily(i,group_index);

            
            casePred = [];
            for k = 1:K_rep
                casePred(:,:,k) = SimNet(sim_span, X0, W_L, psi_t, ZIP, group);
            end
            casePred_mean = mean(casePred,3);   % mean exposure across runs
            for g = 1:4
                group_pos(g,:)  = sum(casePred_mean(group.ZIPs{g},:),1);
            end
            group_pos(5,:)  = sum(casePred_mean(:,:),1);
            group_posRate{group_index}(1:4,i)  = [100*sum(group_pos(1:4,sim_span),2)./mean(group.pop(1:4,sim_span),2)]';
            group_posRate{group_index}(5,i)    = [100*sum(group_pos(5,sim_span),2)/mean(sum(group.pop(:,sim_span),1),2)]';
        end
    end

    figure('Position',[0 0 1000 175])
    tlt = tiledlayout(1, 4);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';
    
    for group_index = 1:4
        nexttile
        for g = 1:5
            for i = 1:length(dhetero)
                posRate{group_index}(g,i) = interp1(hetero{group_index},group_posRate{group_index}(g,:),hetero{group_index}(1)+dhetero(i));
            end
        end
        b = bar([posRate{group_index}(1,:) ; posRate{group_index}(2,:) ; posRate{group_index}(3,:) ; posRate{group_index}(4,:) ; posRate{group_index}(5,:)],'FaceColor','flat','EdgeColor','none');
        for i = 1:length(dhetero)
            for g = 1:5
                b(i).CData(g,:) = hex2rgb(STYLE.colors.dark{g});
            end
        end
        % axis square 
        grid on
        % xlabel('group heterophily (%)')
        ylabel('case rate (%)')
        xlim([0.5 5.5])
        ylim([0 100])
        title([group.name{group_index} ' rewired'])
        xticklabels({'Black' 'Latinx' 'White' 'Mixed' 'City'})
    end

    set(gcf,'renderer','painters');
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize) 
    print(gcf,[parentdir '/figures/structure/groupOpen.eps'],'-depsc')
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
    [parentdir,~,~]=fileparts(pwd);

    K_rep = 10;

    % Without any modifications
    for k = 1:K_rep
        casePred(:,:,k) = SimNet(sim_span, X0, W, psi_t, ZIP, group);
    end
    casePred = mean(casePred,3);
    for g = 1:4
       group_rate(1:4,1,g) = 100*sum(casePred(group.ZIPs{g},:),[1 2])/group.pop_mean(g);
    end
    group_rate(1:4,1,5) = 100*sum(casePred,[1 2])/sum(group.pop_mean);


    for g = 1:4
        index = group.ZIPs{g};
        heterophily.original(1,g) = 100 - 100*sum(W(index,index,:),[1 2 3])/sum(W(index,:,:),[1 2 3]);
    end
    

    for i = 1:4
        % Isolate group i
        W_new = IsolateGroups(ZIP.W,group.ZIPs,[i]);
        for g = 1:4
            index = group.ZIPs{g};
            heterophily.groups(i,g) = 100 - 100*sum(W_new(index,index,:),[1 2 3])/sum(W_new(index,:,:),[1 2 3]);
        end

        for k = 1:K_rep
            casePred(:,:,k) = SimNet(sim_span, X0, W_new, psi_t, ZIP, group);
        end
        casePred = mean(casePred,3);
        for g = 1:4
           group_rate(i,2,g) = 100*sum(casePred(group.ZIPs{g},:),[1 2])/group.pop_mean(g); 
        end
        group_rate(i,2,5) = 100*sum(casePred,[1 2])/sum(group.pop_mean);

        % Isolate nodes in group i
        W_new = IsolateNodes(ZIP.W,group.ZIPs,[i]);
        for g = 1:4
            index = group.ZIPs{g};
            heterophily.nodes(i,g) = 100 - 100*sum(W_new(index,index,:),[1 2 3])/sum(W_new(index,:,:),[1 2 3]);
        end
        for k = 1:K_rep
            casePred(:,:,k) = SimNet(sim_span, X0, W_new, psi_t, ZIP, group);
        end
        casePred = mean(casePred,3);
        for g = 1:4
           group_rate(i,3,g) = 100*sum(casePred(group.ZIPs{g},:),[1 2])/group.pop_mean(g); 
        end   
        group_rate(i,3,5) = 100*sum(casePred,[1 2])/sum(group.pop_mean);
    end    


    figure('Position',[0 0 1000 175])
    tlt = tiledlayout(1, 4);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';

    for i = 1:4
        nexttile
        b = bar(squeeze(group_rate(i,:,:)),'FaceColor','flat','EdgeColor','none');
        for g = 1:4
           b(g).CData = hex2rgb(STYLE.colors.dark{g}); 
        end
        b(5).CData = hex2rgb(STYLE.colors.dark{5});
        % axis square 
        grid on
        title([group.name{i} ' isolated'])
        ylabel('case rate (%)')
        xlim([0.5 3.5])
        ylim([0 100])
        xticklabels({'baseline' 'groups' 'nodes'})
        xtickangle(0)    
    end
    set(gcf,'renderer','painters');
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize) 
    print(gcf,[parentdir '/figures/structure/groupClose.eps'],'-depsc')

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
    print(gcf,[parentdir '/figures/structure/structure.eps'],'-depsc')
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
end













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Trade-offs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Vary beta randomly for each group 
function corr_rb = groupRateBeta(sim_span, X0, ZIP, group, STYLE)
    [parentdir,~,~]=fileparts(pwd);

    figure('Position',[0 0 1000 250])
    tlt = tiledlayout(1, 4);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';



    %% Generate random scalings
    beta_range = [0.05 1];
    psiA_trip_range = [0.5 2];
    N_rand = 100;


    beta_random = rand([4,N_rand])*(beta_range(2)-beta_range(1))+beta_range(1);
    betaScale  = beta_random./repmat(group.beta_mean,1,size(beta_random,2));
    psiScale = betaScale.*rand([4,N_rand])*(psiA_trip_range(2)-psiA_trip_range(1))+psiA_trip_range(1);
    tripScale = sqrt(betaScale./psiScale);

    nexttile()
    hold on
    for i = 1:size(beta_random,2)
        W_hypo    = ZIP.W;
        psi_hypo  = ZIP.psi;
        for g = 1:4
            ZIPs_toRescale = group.ZIPs{g};

            factor_W    = tripScale(g,i);
            factor_psi  = psiScale(g,i);


            W_hypo(ZIPs_toRescale,:,:) = factor_W*W_hypo(ZIPs_toRescale,:,:);
            psi_hypo(ZIPs_toRescale,:) = factor_psi*psi_hypo(ZIPs_toRescale,:);
        end
        hypoStats = HypoStats(sim_span, X0, W_hypo, psi_hypo, ZIP, group);

        x_group = hypoStats.group.beta;
        y_group = hypoStats.group.posRate;
        s_group = group.pop;

        X_group(:,i) = mean(hypoStats.group.beta(:,sim_span),2);
        Y_group(:,i) = hypoStats.group.posRate(:);

    end
    for g = 1:4
        s = scatter(X_group(g,:),Y_group(g,:),40,'filled');
        s.MarkerEdgeColor = 'none';
        s.MarkerFaceColor = STYLE.colors.dark{g};
        s.LineWidth = 1;
    end
    
    index = X_group<0.5;
    xdata = X_group(index);
    ydata = Y_group(index);
    coeff = polyfit(xdata,ydata,1);
    beta_i = linspace(0,1,100);
    plot(beta_i,polyval(coeff,beta_i),'k','LineWidth',2)        
    axis square
    grid on
    xlim([0 1])
    ylim([0 100])
    xlabel('\beta')
    ylabel('case rate (%)')
    
    
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize) 
    set(gcf, 'Renderer', 'painters');
    print(gcf,[parentdir '/figures/tradeoff/illustrateBeta.eps'],'-depsc')
    
    corr_rb.X_group  = X_group;
    corr_rb.Y_group  = Y_group;
    corr_rb.slope    = coeff(1);
    corr_rb.bias     = coeff(2);
end

%% Plot stats/ratio scatters
function Scatter(x_ZIP,y_ZIP,x_group,y_group,group_ZIPs,STYLE)
    hold on
    for g = 1:4
       for i = 1:length(group_ZIPs{g})
           index = group_ZIPs{g}(i);
           scatter(x_ZIP(index),y_ZIP(index),75,hex2rgb(STYLE.colors.light{g}),'filled')
       end
    end
    for g = 1:4
        s = scatter(x_group(g),y_group(g),150,'filled');
        s.MarkerEdgeColor = 'k';
        s.MarkerFaceColor = hex2rgb(STYLE.colors.dark{g});
        s.LineWidth = 1;
    end
    grid on
end




%% Census scatter
function corr_cen = CensusScatter_horizontal(ZIP,group,STYLE)
    [parentdir,~,~]=fileparts(pwd);

    % x_ZIP = [ZIP.census.hhold ZIP.census.income ZIP.census.foreign ZIP.census.ptrans ZIP.census.health ZIP.census.worker ZIP.census.poverty ZIP.census.insur ZIP.census.ocrowd ZIP.census.bo50]';
    % x_group = [group.census.hhold' group.census.income' group.census.foreign' group.census.ptrans' group.census.health' group.census.worker' group.census.poverty' group.census.insur' group.census.ocrowd' group.census.bo50']';
    % scatter_title = {'hhold size' 'income' 'foreign' 'pub trans' 'health' 'worker' 'poverty' 'insurance' 'overcrowded' 'old building'};
    
    
    x_ZIP = [ZIP.census.hhold ZIP.census.income ZIP.census.worker ZIP.census.insur ZIP.census.ocrowd ZIP.census.bo50]';
    y_ZIP = [ZIP.beta_mean; ZIP.posRate];

    x_group = [group.census.hhold' group.census.income' group.census.worker' group.census.insur' group.census.ocrowd' group.census.bo50']';
    y_group = [group.beta_mean group.posRate']';

    scatter_title = {'household size' 'median income' 'employed' 'insured' 'overcrowded' 'old building'};
    fig_title = {'\beta' 'case rate (%)'};
    
    xlabel_text = {'#' '$10k' '%' '%' '%' '%'};

    figure('Position',[0 0 1000 500])
    tlt = tiledlayout(3, 4);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';
        
    for j = 1:length(scatter_title)
       for i = 1:2
           nexttile      
           coeff = polyfit(x_ZIP(j,:),y_ZIP(i,:),1);
           slope(i,j) = coeff(1);
           bias(i,j)  = coeff(2);
           
           xi = linspace(0.9*min(x_ZIP(j,:)),1.1*max(x_ZIP(j,:)),100);
           
           if j == 2
                Scatter(x_ZIP(j,:)/1e4,y_ZIP(i,:),x_group(j,:)/1e4,y_group(i,:),group.ZIPs,STYLE)
                hold on
                plot(xi/1e4,polyval(coeff,xi),'k','LineWidth',2)
           else
                Scatter(x_ZIP(j,:),y_ZIP(i,:),x_group(j,:),y_group(i,:),group.ZIPs,STYLE)
                hold on
                plot(xi,polyval(coeff,xi),'k','LineWidth',2)
           end
        
           [r p] = corrcoef(x_ZIP(j,:)',y_ZIP(i,:)');
           title(scatter_title{j})
           ylabel(fig_title{i})
           
           R(i,j) = r(1,2);
           P(i,j) = p(1,2);
           xlabel(xlabel_text{j})
        end

    end
    
    set(gcf,'renderer','painters');
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize) 
    print(gcf,[parentdir '/figures/tradeoff/censusScatters_horizontal.eps'],'-depsc')
    
    corr_cen.R = R;
    corr_cen.P = P;
    corr_cen.x_ZIP = x_ZIP;
    corr_cen.y_ZIP = y_ZIP;
    corr_cen.x_group = x_group;
    corr_cen.y_group = y_group;
    corr_cen.beta_ZIP = y_ZIP(1,:);
    corr_cen.slope = slope;
    corr_cen.bias = bias;
end


%% How should groups change to achieve a certain outcome
function CensusTripBars(ZIP,group,corr_rb,corr_cen,STYLE)
    [parentdir,~,~]=fileparts(pwd);

    figure('Position',[0 0 1000 250])
    tlt = tiledlayout(1, 4);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';

    pcoeff = [corr_rb.slope corr_rb.bias];
    target_vec = [20 30 40 50 60];
    for i = 1:length(target_vec)
        target_xticklabel{i} = num2str(target_vec(i));
    end
    for g = 1:4
        beta_current = group.beta_mean(g);
        for i = 1:length(target_vec)
            target = target_vec(i);
            beta_target   = (target - pcoeff(2))/pcoeff(1);
            dbeta(g,i)  = (beta_target - beta_current)/beta_current;
            dvuln(g,i)  = 100*dbeta(g,i);
            dtrip(g,i)  = sign(dbeta(g,i))*sqrt(abs(dbeta(g,i))*group.tripRate_mean(g));
        end
    end

    
    nexttile
    b = bar(dvuln','FaceColor','flat','EdgeColor','none');
    for g = 1:4
        b(g).CData = hex2rgb(STYLE.colors.dark{g});
    end
    axis square
    grid on
    title('vulnerability')
    ylabel('%.')
    xlim([0.4 5.6])
    xlabel('case rate (%)')
    xticklabels(target_xticklabel)

    nexttile
    b = bar(dtrip','FaceColor','flat','EdgeColor','none');
    for g = 1:4
        b(g).CData = hex2rgb(STYLE.colors.dark{g});
    end
    axis square
    grid on
    title('trip rate')
    ylabel('#')
    xlim([0.4 5.6])
    xticklabels(target_xticklabel)
    xlabel('case rate (%)')
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize) 
    print(gcf,[parentdir '/figures/tradeoff/targetBars_behavior.eps'],'-depsc')
    
    poly_cen = [corr_cen.slope(1,:); corr_cen.bias(1,:)]';


    for g = 1:4

        beta_current = group.beta_mean(g);
        for i = 1:length(target_vec)
            target = target_vec(i);
            beta_target   = (target - pcoeff(2))/pcoeff(1);
            for j = 1:size(poly_cen,1)
                census_target  = (beta_target - poly_cen(j,2))/poly_cen(j,1);
                census_current = corr_cen.x_group(j,g);
                dcensus{g}(j,i) = census_target - census_current;
                TABLE{j}(g,i) = dcensus{g}(j,i);
            end
        end
    end

    census_text = {'household size' 'median income' 'employed' 'insured' 'overcrowded' 'old building'};
    census_unit = {'#' '$10k' '%.' '%.' '%.' '%.'}; 


    figure('Position',[0 0 1000 500])
    tlt = tiledlayout(2, 3);
    tlt.Padding = 'none';
    tlt.TileSpacing = 'compact';
    
    YTICKS{1} = [-2:0.5:2];
    YTICKS{2} = [-8:2:8];
    YTICKS{3} = [-60:20:60];
    YTICKS{4} = [-40:10:40];
    YTICKS{5} = [-40:10:40];
    YTICKS{6} = [-60:20:60];
    
    YTICKLABELS{1} = {'-2' '' '-1' '' '0' '' '1' '' '2'};
    YTICKLABELS{2} = {'-8' '' '-4' '' '0' '' '4' '' '8'};
    YTICKLABELS{3} = {'-60' '-40' '-20' '0' '20' '40' '60'};
    YTICKLABELS{4} = {'-40' '' '-20' '' '0' '' '20' '' '40'};
    YTICKLABELS{5} = {'-40' '' '-20' '' '0' '' '20' '' '40'};
    YTICKLABELS{6} = {'-60' '-40' '-20' '0' '20' '40' '60'};
    
    for i = 1:6
        nexttile
        if i == 2
            b = bar(TABLE{i}'/1e4,'FaceColor','flat','EdgeColor','none');
        else
            b = bar(TABLE{i}','FaceColor','flat','EdgeColor','none');
        end
        for g = 1:4
            b(g).CData = hex2rgb(STYLE.colors.dark{g});
        end
        axis square
        grid on
        title(census_text{i})
        ylabel(census_unit{i})
        xticklabels(target_xticklabel)
        xlabel('case rate (%)')
        xlim([0.4 5.6])
        yticks(YTICKS{i})
        ylim([min(YTICKS{i}) max(YTICKS{i})])
        yticklabels(YTICKLABELS{i})
    end
    
    set(gcf,'renderer','painters');
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize)
    print(gcf,[parentdir '/figures/tradeoff/targetBars_census.eps'],'-depsc')

end

%% Change required to achieve city average case rate
function PlotTradeOffs(ZIP,group,corr_rb,corr_cen,STYLE)
    [parentdir,~,~]=fileparts(pwd);
    
    caseRate_cityMean = (group.posRate*group.pop_mean)/sum(group.pop_mean);
    case_target = caseRate_cityMean;
    pcoeff = [corr_rb.slope corr_rb.bias];
    beta_target   = (case_target - pcoeff(2))/pcoeff(1);
    
    for g = 1:4
        beta_current = group.beta_mean(g);
    
        ratio_needed(g)    = beta_target/beta_current;
        vuln_target(g)     = ratio_needed(g)*group.vuln_mean(g);
        tripRate_target(g) = sqrt(ratio_needed(g))*group.tripRate_mean(g);
    end
    start(1,:) = group.vuln_mean;
    start(2,:) = group.tripRate_mean;
    
    stop(1,:) = vuln_target;
    stop(2,:) = tripRate_target;
    
    
    start(3,:) = group.census.hhold;
    start(4,:) = group.census.income;
    start(5,:) = group.census.worker;
    start(6,:) = group.census.insur;
    start(7,:) = group.census.ocrowd;
    start(8,:) = group.census.bo50;
    
    for i = 1:6
        slope_actual = corr_cen.slope(2,i);
        bias_actual  = corr_cen.bias(2,i);
    
        case_start(i,:) = start(i+2,:)*slope_actual + bias_actual;
        stop(i+2,:) = start(i+2,:) + (case_target-group.posRate)/slope_actual;
    end
    
    scaling_factors = [1 1 1 1e3 1 1 1 1];
    for i = 1:8
        start(i,:) = start(i,:)/scaling_factors(i);
        stop(i,:)  = stop(i,:)/scaling_factors(i);
    end
    
    xlims = [0 0.02; 6 10 ; 0 4 ; 0 100 ; 0 100 ; 0 100 ; 0 40 ; 0 100];
    titles = {'vulnerability' 'trip rate (#)' 'household size (#)' 'income ($10k)' 'employed (%)' 'insured (%)' 'overcrowded (%)' 'old building (%)'};
    
    figure('Position',[0 0 900 250])
    tlt = tiledlayout(1, 8);
    tlt.Padding = 'tight';
    tlt.TileSpacing = 'compact';
    for i = 1:8
        nexttile(i)
        hold on
        for g = 1:4
            plot([start(i,g) stop(i,g)],[g g],'Color',[1 1 1]*0,'LineWidth',1)
            plot(start(i,g),g,'ko','MarkerFaceColor',STYLE.colors.dark{g},'MarkerSize',15)
            plot(stop(i,g),g,'kp','MarkerFaceColor',STYLE.colors.light{g},'MarkerSize',15)
        end
        grid on
        ylim([0.5 4.5]), yticks([1:4 4.5]), yticklabels({''}),  xlim(xlims(i,:))
        set(gca,'XMinorGrid','on')
        title(titles{i}), % xlabel(units{i})
    end
    set(findall(gcf,'-property','FontSize'),'FontSize',STYLE.fontSize)
        
    % change required
    diff = stop-start;
    rel_diff = 100*((stop-start)./start);

    set(gcf,'renderer','painters');
    print(gcf,[parentdir '/figures/tradeoff/tradeoff_to_average.eps'],'-depsc')
end


%% hex2rgb
function [ rgb ] = hex2rgb(hex,range)
% hex2rgb converts hex color values to rgb arrays on the range 0 to 1. 
% 
% 
% * * * * * * * * * * * * * * * * * * * * 
% SYNTAX:
% rgb = hex2rgb(hex) returns rgb color values in an n x 3 array. Values are
%                    scaled from 0 to 1 by default. 
%                    
% rgb = hex2rgb(hex,256) returns RGB values scaled from 0 to 255. 
% 
% 
% * * * * * * * * * * * * * * * * * * * * 
% EXAMPLES: 
% 
% myrgbvalue = hex2rgb('#334D66')
%    = 0.2000    0.3020    0.4000
% 
% 
% myrgbvalue = hex2rgb('334D66')  % <-the # sign is optional 
%    = 0.2000    0.3020    0.4000
% 
%
% myRGBvalue = hex2rgb('#334D66',256)
%    = 51    77   102
% 
% 
% myhexvalues = ['#334D66';'#8099B3';'#CC9933';'#3333E6'];
% myrgbvalues = hex2rgb(myhexvalues)
%    =   0.2000    0.3020    0.4000
%        0.5020    0.6000    0.7020
%        0.8000    0.6000    0.2000
%        0.2000    0.2000    0.9020
% 
% 
% myhexvalues = ['#334D66';'#8099B3';'#CC9933';'#3333E6'];
% myRGBvalues = hex2rgb(myhexvalues,256)
%    =   51    77   102
%       128   153   179
%       204   153    51
%        51    51   230
% 
% HexValsAsACharacterArray = {'#334D66';'#8099B3';'#CC9933';'#3333E6'}; 
% rgbvals = hex2rgb(HexValsAsACharacterArray)
% 
% * * * * * * * * * * * * * * * * * * * * 
% Chad A. Greene, April 2014
%
% Updated August 2014: Functionality remains exactly the same, but it's a
% little more efficient and more robust. Thanks to Stephen Cobeldick for
% the improvement tips. In this update, the documentation now shows that
% the range may be set to 256. This is more intuitive than the previous
% style, which scaled values from 0 to 255 with range set to 255.  Now you
% can enter 256 or 255 for the range, and the answer will be the same--rgb
% values scaled from 0 to 255. Function now also accepts character arrays
% as input. 
% 
% * * * * * * * * * * * * * * * * * * * * 
% See also rgb2hex, dec2hex, hex2num, and ColorSpec. 
% 
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



function violins = violinplot(data, cats, varargin)
%Violinplots plots violin plots of some data and categories
%   VIOLINPLOT(DATA) plots a violin of a double vector DATA
%
%   VIOLINPLOT(DATAMATRIX) plots violins for each column in
%   DATAMATRIX.
%
%   VIOLINPLOT(DATAMATRIX, CATEGORYNAMES) plots violins for each
%   column in DATAMATRIX and labels them according to the names in the
%   cell-of-strings CATEGORYNAMES.
%
%   In the cases above DATA and DATAMATRIX can be a vector or a matrix,
%   respectively, either as is or wrapped in a cell.
%   To produce violins which have one distribution on one half and another
%   one on the other half, DATA and DATAMATRIX have to be cell arrays
%   with two elements, each containing a vector or a matrix. The number of
%   columns of the two data sets has to be the same.
%
%   VIOLINPLOT(DATA, CATEGORIES) where double vector DATA and vector
%   CATEGORIES are of equal length; plots violins for each category in
%   DATA.
%
%   VIOLINPLOT(TABLE), VIOLINPLOT(STRUCT), VIOLINPLOT(DATASET)
%   plots violins for each column in TABLE, each field in STRUCT, and
%   each variable in DATASET. The violins are labeled according to
%   the table/dataset variable name or the struct field name.
%
%   violins = VIOLINPLOT(...) returns an object array of
%   <a href="matlab:help('Violin')">Violin</a> objects.
%
%   VIOLINPLOT(..., 'PARAM1', val1, 'PARAM2', val2, ...)
%   specifies optional name/value pairs for all violins:
%     'Width'        Width of the violin in axis space.
%                    Defaults to 0.3
%     'Bandwidth'    Bandwidth of the kernel density estimate.
%                    Should be between 10% and 40% of the data range.
%     'ViolinColor'  Fill color of the violin area and data points. Accepts
%                    1x3 color vector or nx3 color vector where n = num
%                    groups. In case of two data sets being compared it can 
%                    be an array of up to two cells containing nx3
%                    matrices.
%                    Defaults to the next default color cycle.
%     'ViolinAlpha'  Transparency of the violin area and data points.
%                    Can be either a single scalar value or an array of
%                    up to two cells containing scalar values.
%                    Defaults to 0.3.
%     'EdgeColor'    Color of the violin area outline.
%                    Defaults to [0.5 0.5 0.5]
%     'BoxColor'     Color of the box, whiskers, and the outlines of
%                    the median point and the notch indicators.
%                    Defaults to [0.5 0.5 0.5]
%     'MedianColor'  Fill color of the median and notch indicators.
%                    Defaults to [1 1 1]
%     'ShowData'     Whether to show data points.
%                    Defaults to true
%     'ShowNotches'  Whether to show notch indicators.
%                    Defaults to false
%     'ShowMean'     Whether to show mean indicator
%                    Defaults to false
%     'ShowBox'      Whether to show the box.
%                    Defaults to true
%     'ShowMedian'   Whether to show the median indicator.
%                    Defaults to true
%     'ShowWhiskers' Whether to show the whiskers
%                    Defaults to true
%     'GroupOrder'   Cell of category names in order to be plotted.
%                    Defaults to alphabetical ordering

% Copyright (c) 2016, Bastian Bechtold
% This code is released under the terms of the BSD 3-clause license

hascategories = exist('cats','var') && not(isempty(cats));

%parse the optional grouporder argument
%if it exists parse the categories order
% but also delete it from the arguments passed to Violin
grouporder = {};
idx=find(strcmp(varargin, 'GroupOrder'));
if ~isempty(idx) && numel(varargin)>idx
    if iscell(varargin{idx+1})
        grouporder = varargin{idx+1};
        varargin(idx:idx+1)=[];
    else
        error('Second argument of ''GroupOrder'' optional arg must be a cell of category names')
    end
end

% check and correct the structure of ViolinColor input
idx=find(strcmp(varargin, 'ViolinColor'));
if ~isempty(idx) && iscell(varargin{idx+1})
    if length(varargin{idx+1}(:))>2
        error('ViolinColor input can be at most a two element cell array');
    end
elseif ~isempty(idx) && isnumeric(varargin{idx+1})
    varargin{idx+1} = varargin(idx+1);
end

% check and correct the structure of ViolinAlpha input
idx=find(strcmp(varargin, 'ViolinAlpha'));
if ~isempty(idx) && iscell(varargin{idx+1})
    if length(varargin{idx+1}(:))>2
        error('ViolinAlpha input can be at most a two element cell array');
    end
elseif ~isempty(idx) && isnumeric(varargin{idx+1})
    varargin{idx+1} = varargin(idx+1);
end

% tabular data
if isa(data, 'dataset') || isstruct(data) || istable(data)
    if isa(data, 'dataset')
        colnames = data.Properties.VarNames;
    elseif istable(data)
        colnames = data.Properties.VariableNames;
    elseif isstruct(data)
        colnames = fieldnames(data);
    end
    catnames = {};
    if isempty(grouporder)
        for n=1:length(colnames)
            if isnumeric(data.(colnames{n}))
                catnames = [catnames colnames{n}]; %#ok<*AGROW>
            end
        end
        catnames = sort(catnames);
    else
        for n=1:length(grouporder)
            if isnumeric(data.(grouporder{n}))
                catnames = [catnames grouporder{n}];
            end
        end
    end
    
    for n=1:length(catnames)
        thisData = data.(catnames{n});
        violins(n) = Violin(thisData, n, varargin{:});
    end
    set(gca, 'XTick', 1:length(catnames), 'XTickLabels', catnames);
    set(gca,'Box','on');
    return
elseif iscell(data) && length(data(:))==2 % cell input
    if not(size(data{1},2)==size(data{2},2))
        error('The two input data matrices have to have the same number of columns');
    end
elseif iscell(data) && length(data(:))>2 % cell input
    error('Up to two datasets can be compared');
elseif isnumeric(data) % numeric input   
    % 1D data, one category for each data point
    if hascategories && numel(data) == numel(cats)    
        if isempty(grouporder)
            cats = categorical(cats);
        else
            cats = categorical(cats, grouporder);
        end
        
        catnames = (unique(cats)); % this ignores categories without any data
        catnames_labels = {};
        for n = 1:length(catnames)
            thisCat = catnames(n);
            catnames_labels{n} = char(thisCat);
            thisData = data(cats == thisCat);
            violins(n) = Violin({thisData}, n, varargin{:});
        end
        set(gca, 'XTick', 1:length(catnames), 'XTickLabels', catnames_labels);
        set(gca,'Box','on');
        return
    else
        data = {data};
    end
end

% 1D data, no categories
if not(hascategories) && isvector(data{1})
    violins = Violin(data, 1, varargin{:});
    set(gca, 'XTick', 1);
% 2D data with or without categories
elseif ismatrix(data{1})
    for n=1:size(data{1}, 2)
        thisData = cellfun(@(x)x(:,n),data,'UniformOutput',false);
        violins(n) = Violin(thisData, n, varargin{:});
    end
    set(gca, 'XTick', 1:size(data{1}, 2));
    if hascategories && length(cats) == size(data{1}, 2)
        set(gca, 'XTickLabels', cats);
    end
end

set(gca,'Box','on');

end
