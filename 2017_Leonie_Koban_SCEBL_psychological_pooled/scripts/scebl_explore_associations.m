% List variables
% ------------------------------------------------------------------------
list_variables(scebl_pooleddata)

ynames = {'ExpectancyEffectonPain  (med temp)'  'SocialInfluenceEffectonPain (med temp)' 'CSEffectonPain  (med temp)'};

covnames = {'Gender (1=M, 2=F)' 'Age'};  % 'Experimenter' - but not fully coded yet

predictornames = scebl_pooleddata.Subj_Level.names(8:27);

% Extract variables of interest
% ------------------------------------------------------------------------

% Conditioning variable
C = get_var(scebl_pooleddata, 'InstructionVersion');

% Outcome effects of interest
clear Y
for i = 1:length(ynames)
    Y(:, i) = get_var(scebl_pooleddata, ynames{i});
end

% Covariates
clear covs
for i = 1:length(covnames)
    covs(:, i) = get_var(scebl_pooleddata, covnames{i});
end
covs = scale(covs, 1); % mean-center

% Predictor effects of interest
clear X
for i = 1:length(predictornames)
    X(:, i) = get_var(scebl_pooleddata, predictornames{i});
end

%% Correlations
% ------------------------------------------------------------------------

[r,t,p,fdrp, fdrthresh] = correlation('r', [Y covs X]);

create_figure('Correlations', 1, 2);
imagesc(r); colorbar
set(gca, 'YDir', 'reverse'); axis tight
title('r, first 3 rows are pain mod fx');

r(isnan(fdrp) | fdrp == 0) = 0;

subplot(1, 2, 2);
imagesc(r, [-1 1]); colorbar
set(gca, 'YDir', 'reverse'); axis tight
title('signif. FDR-corr r, first 3 rows are pain mod fx');

names = [ynames covnames predictornames];
set(gca, 'YTick', 1:length(names), 'YTickLabel', names);

cm = [colormap_tor([.3 0 1], [1 1 1]); colormap_tor([1 1 1], [1 .7 0])];
colormap(cm)

%% Scale data: rank within study

Xr = scale_within_condition(X, C); % normalized rank within study
Yr = scale_within_condition(Y, C);
covsr = scale_within_condition(covs, C, @(X) scale(X, 1)); % mean-center within study

%% Correlations
% ------------------------------------------------------------------------


[r,t,p,fdrp, fdrthresh] = correlation('r', [Yr covsr Xr]);

create_figure('Correlations', 1, 2);
imagesc(r); colorbar
set(gca, 'YDir', 'reverse'); axis tight
title('r, first 3 rows are pain mod fx');

r(isnan(fdrp) | fdrp == 0) = 0;

subplot(1, 2, 2);
imagesc(r, [-1 1]); colorbar
set(gca, 'YDir', 'reverse'); axis tight
title('signif. FDR-corr r, first 3 rows are pain mod fx');

names = [ynames covnames predictornames];
set(gca, 'YTick', 1:length(names), 'YTickLabel', names);

cm = [colormap_tor([.3 0 1], [1 1 1]); colormap_tor([1 1 1], [1 .7 0])];
colormap(cm)

%% Stepwise regression
%
% May be somewhat optimistic, but good as exploratory tool

create_figure('stepwise', 1, length(ynames));

XXnames = [predictornames(1:13) covnames];

for k = 1:length(ynames)
    
    subplot(1, length(ynames), k);
    
    % Must use subset of vars to avoid dropping too many observations
    [wasnan, XX, YY] = nanremove([X(:, 1:13) covs], Y(:, k));
    
    [B,SE,PVAL,INMODEL,STATS,NEXTSTEP,HISTORY] = stepwisefit(XX, YY, 'penter', .1, 'premove', .11);
    
    YYname = format_strings_for_legend(ynames{k});
    YYname = YYname{1};
    
    fprintf('Outcome: %s\nSignificant predictors: ', YYname)
    disp(XXnames(INMODEL))
    
    yhat = XX(:, INMODEL) * B(INMODEL); % does not include intercept
    plot_correlation_samefig(yhat, YY);
    
    set(gca, 'FontSize', 12);
    axis tight
    xlabel(sprintf('Predicted %s', YYname));
    ylabel(YYname);
    
    drawnow
    
end

%% LASSO trace - cross-validated regression

create_figure('LASSO', length(ynames), 3);

for k = 1:length(ynames)
    
    % Must use subset of vars to avoid dropping too many observations
    [wasnan, XX, YY] = nanremove([X(:, 1:13) covs], Y(:, k));
    
    YYname = format_strings_for_legend(ynames{k});
    YYname = YYname{1};
    
    [B,S] = lasso(XX, YY,'CV',10,'PredictorNames', XXnames);
    
    subplot(length(ynames), 3, 3*(k-1) + 1)
    
    axTrace = lassoPlot(B,S, 'PredictorNames', S.PredictorNames, 'Parent', gca);
    %legend(S.PredictorNames)
    % hh = findobj(gca, 'Type', 'Line');
    % set(hh, 'LineWidth', 3);
    
    subplot(length(ynames), 3, 3*(k-1) + 2)
    axTrace = lassoPlot(B,S, 'PredictorNames', S.PredictorNames, 'PlotType', 'CV', 'Parent', gca);
    
    b_inmodel = B(:, S.IndexMinMSE); % cross-validated min MSE
    wh_inmodel = logical(b_inmodel);
    
    subplot(length(ynames), 3, 3*(k-1) + 3)
    
    Xrefit = [XX(:, wh_inmodel) ones(size(XX, 1), 1)];
    bb = pinv(Xrefit) * YY;
    yhat = Xrefit * bb;
    plot_correlation_samefig(yhat, YY);
    
    set(gca, 'FontSize', 12);
    axis tight
    xlabel(sprintf('Predicted %s', YYname));
    ylabel(YYname);
    
    drawnow
    
    fprintf('Outcome: %s\nSignificant predictors: ', YYname)
    disp(XXnames(wh_inmodel))
    
    Coeff = b_inmodel(wh_inmodel);
    Name = XXnames(wh_inmodel)';
    sigtable = table(Name, Coeff);
    disp(sigtable);
    
end
