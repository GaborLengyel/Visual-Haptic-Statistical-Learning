function explicit = Fig4_explicit(Data)

t1 = Data.exp1_V2H.TrialData.explicit;
t2 = Data.exp2_H2V.TrialData.explicit;
data = [t1;t2];
data = [data(~isnan(data{:,3}),:), table([ones(sum(~isnan(t1{:,3})),1); ones(sum(~isnan(t2{:,3})),1)+1], 'VariableNames', {'experiment'})];

data.explicitrel=data.explicit;
data.explicitrel(data.experiment==1)=data.explicitrel(data.experiment==1)/6;
data.explicitrel(data.experiment==2)=data.explicitrel(data.experiment==2)/4;

N=size(data,1);
N1=sum(data.experiment==1);
N2=sum(data.experiment==2);
explicit = zeros(N,4);
%% Controlling for the explicit knowledge

Xr=[ones(N,1) data.haptic data.experiment-1];
Br=regress(data.visual,Xr);

[R_rawj,P_rawj,RLO_rawj,RUP_rawj]=corrcoef([Xr*Br, data.visual]);
fprintf('\n\n  Controlling for explicit knowledge:\n')
fprintf('raw correlation (joint, with separate intercepts): R=%.2g (%.2g-%.2g), p=%.2g\n',R_rawj(2,1),RLO_rawj(2,1),RUP_rawj(2,1),P_rawj(2,1));
explicit(:,1:2) = [Xr*Br, data.visual];

b_rawj=Br(2);
a1_rawj=Br(1);
a2_rawj=Br(1)+Br(3);

ix1=data.experiment==1;
ix2=data.experiment==2;

X1e=[ones(sum(ix1),1) data.explicitrel(ix1)];
X2e=[ones(sum(ix2),1) data.explicitrel(ix2)];

B1v=regress(data.visual(ix1),X1e);
B1h=regress(data.haptic(ix1),X1e);
B2v=regress(data.visual(ix2),X2e);
B2h=regress(data.haptic(ix2),X2e);

[data.visualcorrected3,data.hapticcorrected3]=deal(zeros(size(data.explicitrel)));
data.visualcorrected3(ix1)=data.visual(ix1)-B1v(2)*data.explicitrel(ix1);
data.visualcorrected3(ix2)=data.visual(ix2)-B2v(2)*data.explicitrel(ix2);
data.hapticcorrected3(ix1)=data.haptic(ix1)-B1h(2)*data.explicitrel(ix1);
data.hapticcorrected3(ix2)=data.haptic(ix2)-B2h(2)*data.explicitrel(ix2);

Xe=[ones(N,1) data.hapticcorrected3 data.experiment-1];
Be=regress(data.visualcorrected3,Xe);

[R_corr3j,P_corr3j,RLO_corr3j,RUP_corr3j]=corrcoef([Xe*Be, data.visualcorrected3]);
fprintf('correlation after controlling for explicitness in both training and test (joint, with separate intercepts): R=%.2g (%.2g-%.2g), p=%.2g\n',R_corr3j(2,1),RLO_corr3j(2,1),RUP_corr3j(2,1),P_corr3j(2,1));
explicit(:,3:4) = [Xe*Be, data.visualcorrected3];

b_corr3j=Be(2);
a1_corr3j=Be(1);
a2_corr3j=Be(1)+Be(3);

%% Plotting the distribution of the explicit knowledge

figure(3);
cm=parula(101);
clf;
set(gcf,'Color','white');

cb=axes('Position',[0.06 0.1100 0.9 0.4]);
eu=unique(data.explicitrel);
beu=diag(hist(data.explicitrel,eu));
beu(beu==0)=nan;
hb=bar(cb,eu,beu,'BarWidth',0.7,'BarLayout','Stacked','LineWidth',1);
for i=1:length(eu)
    set(hb(i),'FaceColor',cm(round(eu(i)*100)+1,:),'EdgeColor',[0 0 0],'LineWidth',1,'Clip','on');
end
set(cb,'FontSize',16,'LineWidth',2,'box','off', 'XTick',round(eu*100)/100,'YTick',[0 4 8]);
xlabel('explicitness');
ylabel('number of participants');

%% Plotting the correlation between visual and haptic performance while showing the explicit knowledge of the participants

figure(4);
cm=parula(101);
clf;

set(gcf,'Color','white');

ax=axes('Position',[0.15 0.1100 0.55 0.8150]);
hold on;
plot([0 0],[0.2 1],'-k','LineWidth',2);
plot([-0.4 1],[0.5 0.5],'-k','LineWidth',2);
for i=1:N
    c=cm(round(data.explicitrel(i)*100)+1,:);
    if data.experiment(i)==1
        plot(data.haptic(i),data.visual(i),'o','Color',c,'MarkerFaceColor',c,'MarkerSize',10,'LineWidth',2);
    elseif data.experiment(i)==2
        plot(data.haptic(i),data.visual(i),'s','Color',c,'MarkerFaceColor',c,'MarkerSize',10,'LineWidth',2);
    end
end
lx=[-0.4 (1-a1_rawj)/b_rawj];
plot(lx,a1_rawj+b_rawj*lx,'-r','LineWidth',2);
plot(lx,a2_rawj+b_rawj*lx,':r','LineWidth',2);
hold off;
xlim([-0.4 1]);
ylim([0.2 1]);
set(gca,'FontSize',16,'LineWidth',2,'box','off');
xlabel('haptic pulling performance (\rho)');
ylabel({'visual familiarity performance','(fraction correct)'});
title('raw performance');
%set(gca,'ColorMap',cm);

%% plotting the correlation between visual and haptic performance while controlling for explicit knowledge

figure(5);
cm=parula(101);
clf;
set(gcf,'Color','white');
hold on;
plot([0 0],[0.2 1],'-k','LineWidth',2);
plot([-0.4 1],[0.5 0.5],'-k','LineWidth',2);
for i=1:N
    c=cm(round(data.explicitrel(i)*100)+1,:);
    if data.experiment(i)==1
        plot(data.hapticcorrected3(i),data.visualcorrected3(i),'o','Color','black','MarkerFaceColor','black','MarkerSize',10,'LineWidth',2);
    elseif data.experiment(i)==2
        plot(data.hapticcorrected3(i),data.visualcorrected3(i),'s','Color','black','MarkerFaceColor','black','MarkerSize',10,'LineWidth',2);
    end
end
lx=[-0.4 (1-a1_corr3j)/b_corr3j];
plot(lx,a1_corr3j+b_corr3j*lx,'-r','LineWidth',2);
plot(lx,a2_corr3j+b_corr3j*lx,':r','LineWidth',2);
hold off;
xlim([-0.4 1]);
ylim([0.2 1]);
set(gca,'FontSize',16,'LineWidth',2,'box','off');
xlabel('haptic pulling performance (\rho)');
ylabel({'visual familiarity performance','(fraction correct)'});
title({'residual performance','controlling for explicitness'});
set(gca,'Position',get(ax,'Position'));

end
