function Fig3_forces(Data)

figure(2)
exp_name = fieldnames(Data)';

for experiment=1:length(exp_name)
    
    D = Data.(exp_name{experiment}).forces;
    
    clear ffinal ff
    
    N= D.Trials; %total number of trials across al subjects
    %  ntrial=length(trial); %number of trials for each subject
    %  N=ntrial; %total number of trials across al subjects
    
    ind1 = myfindfirst(D.State'>=9);
    ind2 = min( ind1+3000,myfindfirst(D.State'>=10));
    ind1(ind1==0)=1; %hack for missing data
    ind2(ind2==0)=1;
    
    m=max(ind2-ind1);
    
    pad=m-(ind2-ind1);
    
    %chop data
    chopnan = @(x,ind1,ind2,N) cell2mat(arrayfun(@(k)(cat(4,x(k,1:2,1:2,ind1(k):ind2(k)),NaN(1,2,2,pad(k)))), 1:N, 'UniformOutput', false)');
    
    force    =  chopnan(D.RobotForces,ind1,ind2,N);
    
    gen=strcmp(D.TrialData.blockname,'Generalization');
    s1=gen & D.TrialData.pullvert;
    s2=gen & ~D.TrialData.pullvert;  
    
    bf=D.TrialData.breakforce;
    ub=unique(bf);
    col='rgb';
    
    map=load('rho_map.m');
    
    red=[ 238 34 12]/255;
    green=[ 29 177 0]/255;
    blue=[ 0 118 186]/255;
    
    set(gcf,'Position', [100 100  1200 1000])
     
    subplot(2,2,1+(experiment-1)*2)
    if 0
        q=force(:,1,1,500)>30;
        force(q,1,1,:)=NaN;
        
        q=force(:,1,2,500)>30;
        force(q,1,2,:)=NaN;
    end 
    
    t=linspace(0,3,3001);
    for j=1:length(ub)
        
        for k=1:max(D.Subj)
            s=D.Subj==k;
            ff(k,:)=0.5*(nanmean(squeeze(force(s&s2&bf==ub(j),1,1,:)))+nanmean(squeeze(force(s&s1&bf==ub(j),1,2,:))));
            ffinal(k,j)=0.5*nanmean(D.TrialData.PullForce(s  & gen & D.TrialData.breakforce==ub(j),:));
        end
        
        lw=3;
        if j==1
            [h,g]=shadeplot(t,nanmean(ff),stderr(ff),'r-',red);
            set(h,'Color',red)
            plot(t,(t==t)*7.5,'r--','Color',red,'LineWidth',3)
        elseif (j==2 & experiment==1) |(j==3 & experiment==2)
            [h,g]=shadeplot(t,nanmean(ff),stderr(ff),'b-',blue);
            set(h,'Color',blue)
            plot(t,(t==t)*22.5,'b--','Color',blue,'LineWidth',3)
        elseif experiment==2 & j==2
            [h,g]=shadeplot(t,nanmean(ff),stderr(ff),'b-',green);
            set(h,'Color',green)
            plot(t,(t==t)*15,'g--','Color',green,'LineWidth',3)
            
        end
        set(g,'FaceAlpha',0.3)
        h.LineWidth=lw;
        
    end
       
    xlabel('time (s)')
    ylabel('pulling force (N)')
    axis([0 3 5 36])
    
    shg
    
    subplot(2,2,2+(experiment-1)*2)
    
    colormap(map)
    c=colormap;
    
    for k=1:rows(ffinal)
        s=D.Subj==k;
        f1=0.5*D.TrialData.PullForce( s&gen ,:);
        f2=0.5*D.TrialData.breakforce(s&gen ,:);
        
        cc=corrcoef(f1,f2);
        cc=cc(1,2);
        r=round(32+32*cc);
        plot(ub/2,ffinal(k,:),'Color',c(r,:),'LineWidth',3)
        hold on
    end
    
    plot([5 25],[5 25],'k:','LineWidth',3)
    xlabel('breakage force (N)')
    ylabel('pulling force (N)')
    set(gca,'Clim',[-1 1])
    
    h=colorbar;
    h.Label.String = 'correlation (\rho)';
    set(h,'Limit',[-0.3 1]')
    axis equal
    axis([5 25 5 36])
    nudge_plot(gca,-0.2,0)
    box off
end

text_outside(0.05,0.97,'A. pulling test after visual statistical exposure','FontSize',26,'FontWeight','bold')
text_outside(0.05,0.47,'B. pulling test after haptic statistical exposure','FontSize',26,'FontWeight','bold')

%export_fig FigS3.png
shg

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(gcf,'Fig3.pdf','-dpdf','-r300')
end