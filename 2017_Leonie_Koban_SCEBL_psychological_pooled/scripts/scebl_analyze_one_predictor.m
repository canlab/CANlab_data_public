function scebl_analyze_one_predictor(scebl_pooleddata, xname, yname)

X = get_var(scebl_pooleddata, xname);
Y = get_var(scebl_pooleddata, yname);

condition_on = 'InstructionVersion';

C = get_var(scebl_pooleddata, condition_on);

u = unique(C)';
clear Xi Yi wasnan N Within_study_p

[ri, rirank, Within_study_p, Within_study_rank_p] = deal(NaN * ones(size(u')));
Expt_version = u';

for i = u
   
    wh = C == i;
    
    Xi{i} = X(wh, :);
    Yi{i} = Y(wh, :);
   
    [wasnan{i}, Xi{i}, Yi{i}] = nanremove(Xi{i}, Yi{i});
    
    if ~isempty(Xi{i})
    [ri(i, 1), Within_study_p(i, 1)] = corr(Xi{i}, Yi{i});
    end
    
    N(i, 1) = length(Xi{i});
    
    % "Nonparametric" : Normalized rank data within-study
    % Divide by length (N observations) to keep range same for all studies
    Xirank{i} = rankdata(Xi{i}) ./ length(Xi{i});
    Yirank{i} = rankdata(Yi{i}) ./ length(Yi{i});
   
    if ~isempty(Xi{i})
    [rirank(i, 1), Within_study_rank_p(i, 1)] = corr(Xirank{i}, Yirank{i});
    end
    
end

Within_study_r = ri;
Within_study_rank_r = rirank;

results_table = table(Expt_version, N, Within_study_r, Within_study_p, Within_study_rank_r, Within_study_rank_p);

fprintf('VARIABLES: %s predicting %s\n', xname, yname);
disp(results_table);

%% Version plots

condcolors = seaborn_colors(length(u));

create_figure(xname, 2, 3);

plot_correlation_samefig(X, Y);

for i = 1:length(u)
    plot(Xi{i}, Yi{i}, 'o', 'Color', condcolors{i} ./ 2, 'MarkerFaceColor', condcolors{i});
end
xlabel(xname);
ylabel(yname);
%title('Overall correlation');

subplot(2, 3, 2);

barplot_columns(ri', 'colors', condcolors, 'nofig', 'noviolin', 'noind');
xlabel('Experiment version');
ylabel('Within-study correlation');
title('Within-study correlations');

subplot(2, 3, 3);

barplot_columns(Yi, 'colors', condcolors, 'nofig', 'noviolin');
xlabel('Experiment version');
ylabel(yname);
title('Within-study effects');


subplot(2, 3, 4);

Xrank = cat(1, Xirank{:});
Yrank = cat(1, Yirank{:});

plot_correlation_samefig(Xrank, Yrank);

for i = 1:length(u)
    plot(Xirank{i}, Yirank{i}, 'o', 'Color', condcolors{i} ./ 2, 'MarkerFaceColor', condcolors{i});
end
xlabel(['Normed Rank ' xname]);
ylabel(yname);
%title('Overall within-study rank correlation');

subplot(2, 3, 5);

barplot_columns(rirank', 'colors', condcolors, 'nofig', 'noviolin', 'noind');
xlabel('Experiment version');
ylabel('Within-study correlation');
title('Within-study rank correlation');

subplot(2, 3, 6);
axis off
% barplot_columns(Yirank, 'colors', condcolors, 'nofig', 'noviolin');
% xlabel('Experiment version');
% ylabel(yname);
% title('Within-study effects');


end % function


