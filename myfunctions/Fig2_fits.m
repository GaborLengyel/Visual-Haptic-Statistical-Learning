function results = Fig2_fits(Data)

figure(1)
set(gcf,'Position', [100 100  1200 500])
exp_name = fieldnames(Data)';

for experiment=exp_name
    
    if strcmp(experiment{1}, exp_name{2})
        subplot(1,2,2);
    else
        subplot(1,2,1);
    end
    clear pval psl m fam slc slh slv hcorr O
    
    TrialData = Data.(experiment{1}).TrialData;
    haptics.TrialData=Data.(experiment{1}).TrialData.haptic;
    
    %find clipping
    nsubj=max(haptics.TrialData.Subj);
    for k=1:nsubj
        s=haptics.TrialData.Subj==k & haptics.TrialData.Clamp==1;
        
        f=haptics.TrialData.PullForce(s);
        m(k)=mean(f>79.9);
    end
    fprintf('Percent of trials which saturate %2.2f%%\n',100*mean(m))
    
    %%  caculate the visual familiarity (in fraction correct) and the haptic pulling (correlation) performance 
    for k=1:nsubj
        s=TrialData.familiarity.Subj==k;
        tmp2=TrialData.familiarity(s,:);
        corr=(tmp2.truepair1 & tmp2.Choice==1) |   (tmp2.truepair2 & tmp2.Choice==2);
        fam(k)=mean(corr);
        numtrial=rows(tmp2);
        
        if strcmp(experiment{1}, exp_name{2})
            s=TrialData.haptic.Subj==k  & TrialData.haptic.block_index==2;
        else
            s=TrialData.haptic.Subj==k ;
        end
        
        tmp=haptics.TrialData(s,:);
        t=polyfit(tmp.breakforce,tmp.PullForce,1);
        slh(k)=t(1);
        t=corrcoef(tmp.breakforce,tmp.PullForce);
        hcorr(k)=t(1,2);
        name{k}=tmp.subjname{1}(1:3);
    end
    
    x=hcorr';
    y=fam';
    vis.(experiment{1}) = y;
    hap.(experiment{1}) = x;
    O=table;
    O.name=name';
    O.hcorr=hcorr';
    O.fam=fam';
    
    %% Find the optimal fit for the rectified exponential-binomial model
    N=length(x);
    
    Ru=myweibull_fit.weibull_fit(x,y,numtrial,[-100 -100 0],[100 100 100]); % the function class for the rectified exponential-binomial curve fitting
    
    Y0=mean(y)*ones(N,1);
    
    opt_logl=numtrial*sum(cross_surprise(y,Ru.yfit)+cross_surprise(1-y,1-Ru.yfit)); % the likelihood for the rectified exponential-binomial model
    red_logl=numtrial*sum(cross_surprise(y,Y0)+cross_surprise(1-y,1-Y0)); % the likelihood for the null model  (straight line at the mean familiarity level)
    
    % likelihood ratio test
    chi=2*(opt_logl-red_logl);
    p=1-chi2cdf(chi,2);
    err = sqrt(diag(inv(Ru.hess)));
    
    %log bayes factors for comparing the two models
    b_full= Ru.logl + log(2*pi)*length(Ru.bhat)/2 -0.5*logdet(Ru.hess);
    
    syms p0 ys real
    
    nlogl0 = -( ys * log(p0) +  (1-ys) * log(1-p0));
   
    H0=hessian(nlogl0,  p0);
    f=matlabFunction(H0);
    h0=numtrial*sum(f(Y0,y));
    b_red=red_logl + log(2*pi)*1/2 -0.5*logdet(h0);
    bf=b_full-b_red;
   
    
    %% bootstrap to get 95%CI arounf the maximum likelihood fist using the profile likelihood method
    clear nnl
    N=10^5;
    samp=mvnrnd(repmat(Ru.bhat,N,1),inv(Ru.hess));
    
    for k=1:N
        nnl(k)=myweibull_fit.weibull_model(samp(k,:),x,y,numtrial,0);
    end
    
    dlogl=nnl+Ru.logl;
    s=dlogl <chi2inv(0.95,3)/2;
    
    samp=samp(s,:);
    p0=samp(:,1);
    p1=samp(:,2);
    lambda=samp(:,3);
    
    s=Ru.xi>0;
    ypred=zeros(rows(samp),length(Ru.xi));
    
    p0=1./(1+exp(p0));
    p1=1./(1+exp(p1));
    
    tmp=arrayfun(@(x) p1+(p0-p1).*(exp(-(x./lambda))),Ru.xi(s), 'UniformOutput', false);
    ypred(:,s)=cell2mat(tmp);
    
    tmp=arrayfun(@(x) p0,Ru.xi(~s), 'UniformOutput', false);
    ypred(:,~s)=cell2mat(tmp);
    
    uypredp=[ypred;Ru.yi];
    
    s=sort(ypred);
    lo=s(1,:);
    hi=s(end,:);
    %% Plotting Figure 2
    red=[ 238 34 12]/255;
    green=[ 29 177 0]/255;
    blue=[ 0 118 186]/255;
    set(0,'DefaultAxesFontSize',20);
    map=load('rho_map.m');
    
    colormap(map);
    c=colormap;
    hold on
    
    plot(Ru.xi,Ru.yi,'r','LineWidth',2,'Color',red);
    
    h=fill([Ru.xi fliplr(Ru.xi)],[lo fliplr(hi)],red);
    set(h,'FaceAlpha',.3,'LineStyle','none');
    
    hold on
    set(gca,'clim',[-1 1])
    dx = 0.01; dy = 0.0; % displacement so the text does not overlay the data points
    
    if strcmp(experiment{1}, exp_name{1})
        colorbar off
    end
      
    axis([-0.8 1.01 0.25 1.0]);
    aa=axis;
    plot(aa(1:2),[0.5 0.5],'k-');
    plot([0 0],aa(3:4),'k-'); 
    
    plot(0,mean(y),'go','MarkerFaceColor',green,'MarkerEdgeColor',green,'MarkerSize',15);
    errbar(0,mean(y),1.96*stderr(y),'g+-','LineWidth',5,'Color',green);
    plot(mean(x),0.5,'bo','MarkerFaceColor',blue,'MarkerEdgeColor',blue,'MarkerSize',15);
    errbar(mean(x),0.5,1.96*stderr(x),'horiz','b+-','LineWidth',5,'Color',blue);
    set(gca,'Ytick',0.25:0.25:1);
    set(gca,'Xtick',-0.75:0.25:1);
    scatter(x,y,50,'k','filled')
    
    box off;
    xlabel('haptic pulling performance (\rho)');
    if strcmp(experiment{1}, exp_name{1})
        ylabel({'visual familiarity performance','(fraction correct)'});
    end
    shg
    
    %
    [h,px,~,statx]=ttest(x);
    [h,py,~,staty]=ttest(y-0.5);
    
    
    % Printing out the results
    fprintf('\n\n  Main results for %s:\n', experiment{1})
    fprintf('Haptic pulling performance: Rho= %2.2f %2.2f  with t(%i)=%2.2f  p=%e\n', mean(x),stderr(x),statx.df,statx.tstat,px)
    fprintf('Visual familiarity performance: fraction correct= %2.2f %2.2f  with t(%i)=%2.2f p=%e\n', mean(y),stderr(y),staty.df,staty.tstat,py)
    fprintf('Rectified exponential-binomial fit: likelihood ratio test p=%e chi(2)=%2.2f, log Bayes factor %2.2f, Range convered = %2.2f%%\n',p,chi,bf,100*range(Ru.yi)/0.5)
    
    %%
    % Analysing Within-subject consistency
    cons_corr = zeros(1,nsubj);
    for k=1:nsubj
        
        s=TrialData.familiarity.Subj==k;
        tmp2=TrialData.familiarity(s,:);
        % getting the true pairs out from the familiarity trials
        true_pairs_trials = tmp2.symbols1(logical(tmp2.truepair1),:);
        true_pairs = unique(true_pairs_trials,'rows');

        if strcmp(experiment{1}, exp_name{2})
            s=TrialData.haptic.Subj==k  & TrialData.haptic.block_index==2;
        else
            s=TrialData.haptic.Subj==k ;
        end

        tmp=haptics.TrialData(s,:);
        
        npairs = size(true_pairs,1);
        
        visual_score = zeros(1,npairs);
        haptic_score = zeros(1,npairs);
        for tp=1:npairs

            tmp2_pair = tmp2((((tmp2.symbols1(:,1)==true_pairs(tp,1))+(tmp2.symbols1(:,2)==true_pairs(tp,2)))==2) | (((tmp2.symbols2(:,1)==true_pairs(tp,1)) + (tmp2.symbols2(:,2)==true_pairs(tp,2)))==2),:);
            visual_score(tp)= mean((tmp2_pair.truepair1 & tmp2_pair.Choice==1) | (tmp2_pair.truepair2 & tmp2_pair.Choice==2));

            tmp_pair = tmp((((tmp.symbols1(:,1)==true_pairs(tp,1))+(tmp.symbols1(:,2)==true_pairs(tp,2)))==2) | (((tmp.symbols2(:,1)==true_pairs(tp,1)) + (tmp.symbols2(:,2)==true_pairs(tp,2)))==2),:);
            haptic_score(tp)= mean(tmp_pair.PullForce(tmp_pair.breakforce==45)) - mean(tmp_pair.PullForce((tmp_pair.breakforce==15)&(tmp_pair.truepair1==1)&(tmp_pair.truepair2==1)));

        end

        t=corrcoef(visual_score,haptic_score);
        cons_corr(k)=t(1,2);

    end

cons.(experiment{1}) = cons_corr;

end

% printing out the result of the within-subject consistency analysis
c = [cons.(exp_name{1}), cons.(exp_name{2})];
c = c(~isnan(c));
[H,P,~,statx] = ttest(c);
fprintf('\n\n  Within-subject consistency:\n')
fprintf('Correlation between familiarity and pulling force for individual scenes: Rho= %2.2f %2.2f  with t(%i)=%2.2f  p=%e\n', mean(c),stderr(c),statx.df,statx.tstat,P)

%% Comparing the two experiments
fprintf('\n\n  Comparing the two exeriments:\n')
[H,P,~,statx] = ttest2(vis.(exp_name{1}),vis.(exp_name{2})); % comparing the visual familiarity performance across experiments
fprintf('Difference between the visual familiarity performance in the two experiment: t(%i)= %2.2f,  p=%e\n', statx.df,statx.tstat,P)
[H,P,~,statx] = ttest2(hap.(exp_name{1}),hap.(exp_name{2})); % comparing the haptic pulling performance across experiments
fprintf('Difference between the haptic performance in the two experiment: t(%i)= %2.2f,  p=%e\n', statx.df,statx.tstat,P)
results = [vis.(exp_name{1}), vis.(exp_name{2}), hap.(exp_name{1}), hap.(exp_name{2}), cons.(exp_name{1})', cons.(exp_name{2})'];

%%
shg

%% Plotting Figure 2 cont.

for k=1:2
    subplot(1,2,k)
    nudge_plot(gca,0,0.02,1,0.90)
end
text_outside(0.2,0.97,'visual exposure','FontSize',26,'FontWeight','bold');
text_outside(0.65,0.97,'haptic exposure','FontSize',26,'FontWeight','bold');

text_outside(0.04,0.97,'A','FontSize',26,'FontWeight','bold');
text_outside(0.53,0.97,'B','FontSize',26,'FontWeight','bold');

%export_fig  fig2.png

% print(gcf, '-dpdf', 'Fig2.pdf');

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

%print(gcf,'Fig2.pdf','-dpdf','-r0')

end
