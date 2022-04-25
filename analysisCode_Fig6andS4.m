%% load file names 
clear
case1=["1","2","3","4","5"]; %case number for different conditions
frac=["0","000003125","00000625","0000125","nodelay","00005","0001","0002"];
%naming for different fraction(s) of fast elongating F-actin(%)

for ii=1:5
for jj=5:8

%% for motor contraction, for Figure 6F
data2=load(sprintf('MotCenPos_%s_%s',frac(jj),case1(ii)));
data2(1:6,:)=[]; %data2(1:6,:)=[]; for delays

[col,row]=find(data2(:,5)==0);
timestepm=numel(col);
col(timestepm+1)=numel(data2)/5+1;
for t=1:timestepm
    ini=col(t)+1;final=col(t+1)-1;
    Nmotor(t)=final-ini;
    
x_data=data2(ini:final,2)*0.14; %read x position
y_data=data2(ini:final,3)*0.14; %read y coordinates
x_data(x_data>20)=x_data(x_data>20)-20;x_data(x_data<0)=x_data(x_data<0)+20; %PBC
y_data(y_data>20)=y_data(y_data>20)-20;y_data(y_data<0)=y_data(y_data<0)+20; %PBC

xcenter=mean(x_data);ycenter=mean(y_data);
r=sqrt((x_data-xcenter).^2+(y_data-ycenter).^2);
radius(t)=mean(r);

figure(1)
radius(radius==0)=NaN;
contraction=(radius-radius(1))/radius(1);
contraction(contraction>0)=NaN;
end

%% Figure S4D
plot((1:timestepm),abs(contraction(1:timestepm)),'LineWidth',2)
a=gca;
a.FontSize=20;
a.LineWidth=2;
xlabel('Time [s]')
ylabel('Motor contraction')
box on
saveas(gcf,sprintf('Motorcontract2_%s%s.tiff',frac(jj),case1(ii)))
motorc(ii,jj)=max(abs(contraction(1:15))); %data for motor contraction (maximal)
close
clear data2 col radius contraction

end
end

%% plot motor contraction (Fraction), Figure S4D
y=mean(motorc);z=std(motorc);
figure(2)
errorbar(1:8,y,z,'b-^','LineWidth',2,'MarkerSize',15')
a=gca;
a.FontSize=20;
a.LineWidth=2;
a.XTick=[1:8];
a.XTickLabel={'0','0.03','0.06','0.13','0.25','0.5','1','2'};
xlabel('Fraction of fast elongating F-actin [%] ')
ylabel('Motor contraction')
xlim([1 8])
box off

%% Calculating actin contraction, for Figures 6B,E
clear
case1=["1","2","3","4","5"];
frac=["0","000003125","00000625","0000125","000025","00005","0001","0002"];

for ii=1:5
for jj=1:8

%% for motor contraction
data2=load(sprintf('MotCenPos_f%s_%s',frac(jj),case1(ii)));
data2(1:6,:)=[]; %data2(1:6,:)=[]; for delays

[col,row]=find(data2(:,5)==0);
timestepm=numel(col);
col(timestepm+1)=numel(data2)/5+1;
t=1;
    ini=col(t)+1;final=col(t+1)-1;
    Nmotor(t)=final-ini;
    
x_data=data2(ini:final,2)*0.14; %read x position
y_data=data2(ini:final,3)*0.14; %read y coordinates
a1=min(x_data);a2=max(x_data);a3=min(y_data);a4=max(y_data);

a1=a1-1;a2=a2+1;a3=a3-1;a4=a4+1;
clear x_data y_data timstepm t row radius Nmotor ini final data2 contraction col a
data=load(sprintf('MlbActPos_f%s_%s',frac(jj),case1(ii))); %load actin positions
Nactin=data(1,1);

timestep=(numel(data(:,1))-1)/Nactin;
data1=data;data1(1,:)=[];
nbins_cluster=20;
for t=1:timestep
ini=1+Nactin*(t-1);
final=Nactin*t;

x_data=data1(ini:final,2)*0.14; %read x position
y_data=data1(ini:final,3)*0.14; %read y coordinates
x_data(x_data>20)=x_data(x_data>20)-20;x_data(x_data<0)=x_data(x_data<0)+20; %PBC
y_data(y_data>20)=y_data(y_data>20)-20;y_data(y_data<0)=y_data(y_data<0)+20; %PBC
pos=find(x_data>a1&x_data<a2&y_data>a3&y_data<a4);
x_data=x_data(pos);y_data=y_data(pos);
 
h = histogram2(x_data,y_data,nbins_cluster);
    intensity=h.Values;
intensity=reshape(intensity, [1 numel(intensity)]);
meaninten(t)=mean(intensity);
hetero(t)=std(intensity)/meaninten(1);

clear x_data y_data  intensity std_col std_row xdis ydis dis pos
end
close
figure(1)
plot(1:timestep,hetero(1)./hetero,'LineWidth',2)
a=gca;
a.FontSize=20;
a.LineWidth=2;
xlabel('Time [s]')
ylabel('Actin contraction')
box on
set(gcf,'Position',[1145 523 603 442])

saveas(gcf,sprintf('Actincontract_d%s_%s.tiff',frac(jj),case1(ii)));
close

result(1)=1-min(hetero(1)./hetero);
result(2)=1-hetero(1)/mean(hetero(30:32));
result(3)=max(hetero-hetero(1));
result(4)=(mean(hetero(timestep)-hetero(1)));

maxactin(ii,jj)=result(1);
meanactin(ii,jj)=result(2);
meanactin(meanactin<0)=0;
maxactinsubtract(ii,jj)=result(3);
meanactinsub(ii,jj)=result(4);

clear hetero data1 data2 result
end
end


%% plot actin contraction, Figure 6B, for Figure 6C, replace data with force data
figure(1)
plot(1:8,maxactin,'b-^','LineWidth',2,'MarkerSize',15')
hold on
plot(1:8,abs(meanactin),'r-o','LineWidth',2,'MarkerSize',15')
set(gcf,'Position',[1145 523 603 442])

a=gca;
a.FontSize=20;
a.LineWidth=2;
a.XTick=[1:8];
a.XTickLabel={'0','0.03','0.06','0.13','0.25','0.5','1','2'};
xlabel('Fraction of fast elongating F-actin [%] ')
ylabel('Actin contraction')
ylim([0 1])
xlim([1 8])
legend('Maximum','Plateau')
legend boxoff
box off


%% find force
clear
% file naming and loading
case1=["1","2","3","4","5"];
%frac=["0","000003125","00000625","0000125","000025","00005","0001","0002"];
frac=["0","000003125","00000625","0000125","000025","ft5fa10","nodelay","0002"];

for ii=1:5
for jj=1:8 

data3=load(sprintf('InfoIndv_f%s_%s',frac(jj),case1(ii)));
Nactin=data3(1,1);
timestep=(numel(data3(:,1))-1)/Nactin;
data3(1,:)=[];

for t=1:timestep
    
ini=1+Nactin*(t-1);
final=Nactin*t;
x_data=data3(ini:final,5); %read x position
y_data=data3(ini:final,6); %read y coordinates
numberm(t)=nnz(find(y_data>-1));
%if t>10
%    y_data=data3((1+Nactin*9):(Nactin*10),6);
%end
x_data=data3(ini:final,5); %read x position
fast=x_data(find(y_data>-1));
nonfast=x_data(find(y_data==-1));
fast2=nonfast;fast2(1:numel(fast2))=NaN;
if t==10
    pos=find(y_data>-1);
  
end
fast2(1:numel(fast))=fast;
if t<10
force(t)=sum(abs(fast))/1000;forceall(t)=(sum(abs(fast))+sum(abs(nonfast)))/1000;
end
if t>=10
      fast1=x_data(pos);
 force(t)=sum(abs(fast1))/1000;forceall(t)=(sum(abs(fast1))+sum(abs(nonfast)))/1000;
end
fast0=fast;
clear fast nonfast fast2 x_data y_data fast1
end
plot(0:timestep-1,force/111,'LineWidth',2)
xlabel('Time [s]')
ylabel('Sum of force [nN]')
box on
a=gca;
a.FontSize=18;
a.LineWidth=2;

saveas(gcf,sprintf('sumforce_%s.tiff',frac(jj),case1(ii)));
close

forcesum(ii,jj)=max(force);
numberfast(ii,jj)=max(forceall);
clear data3 force

end
end
%% end of finding sum of tensile forces


%% plot pointed end, Figures 6E-G, use precalculated dataset(s). 
figure(1)
%force0=[0.924217817	1.73618509];
bar(1:2,force0(1:2),'BarWidth',0.5)
a=gca;
a.FontSize=20;
a.LineWidth=2;
ylabel('Actin contraction')
ylim([0 2.0])
set(gcf,'Position',[766,299,328,566])
hold on
errorbar(1:2,[0.924217817	1.73618509],[0.067148072	0.17971013],'r.')
set(a,'FontSize',20,'LineWidth',2,'XTick',[1 2],'XTickLabel',...
    {'Control','Pointed end'},'XTickLabelRotation',45);
box off

figure(2)
%force0=[0.313707348	0.585880234];
bar(1:2,force0(1:2),'BarWidth',0.5)
a=gca;
a.FontSize=20;
a.LineWidth=2;
ylabel('Motor contraction')
ylim([0 1])
hold on
errorbar(1:2,[0.313707348	0.585880234],[0.07175231	0.05380344],'r.')
set(gcf,'Position',[766,299,328,566])
set(a,'FontSize',20,'LineWidth',2,'XTick',[1 2],'XTickLabel',...
    {'Control','Pointed end'},'XTickLabelRotation',45);
box off

figure3=figure;
axes1 = axes('Parent',figure3);

%force0=[80.96047845	41.51199574];
factor1=[106.7 111];
bar(1:2,force0(1:2)./factor1,'BarWidth',0.5)
a=gca;
a.FontSize=20;
a.LineWidth=2;
ylabel('Sum of forces [nN]')
ylim([0 1])
hold on
errorbar(1:2,[80.96047845	41.51199574]./factor1,[5.159350958	1.85671283]./factor1,'r.')
%set(gcf,'Position',[766,445,328,420])
set(gcf,'Position',[766,299,328,566])

set(axes1,'FontSize',20,'LineWidth',2,'XTick',[1 2],'XTickLabel',...
    {'Control','Pointed end'},'XTickLabelRotation',45);
box off

figure3=figure;
axes1 = axes('Parent',figure3);

%force0=[91.4745 47.9438];
factor1=[106.7 111];
bar(1:2,force0(1:2)./factor1,'BarWidth',0.5)
a=gca;
a.FontSize=20;
a.LineWidth=2;
ylabel('Sum of forces [nN]')
ylim([0 1])

%set(gcf,'Position',[766,445,328,420])
set(gcf,'Position',[766,299,328,566])

set(axes1,'FontSize',20,'LineWidth',2,'XTick',[1 2],'XTickLabel',...
    {'Control','Pointed end'},'XTickLabelRotation',45);

%% end of plots