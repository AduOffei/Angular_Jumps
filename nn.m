dp =importdata("../nf100_100/gd_angle_vs_time_dp_1.mat");
hh =importdata("../nf100_100/gd_angle_vs_time_hh_1.mat");
pdbstruct = pdbread('fluc.pdb');
%dp = dp(dp(:,1)>60,:);
%hh = hh(hh(:,1)>60,:);
yel = [dp(dp(:,1)>60,:) ;hh(hh(:,1)>60,:)];
yel40 = [dp(dp(:,1)>40,:) ;hh(hh(:,1)>40,:)];
yello = [dp(dp(:,1)<20,:); hh(hh(:,1)<20,:)];


%%%%%%%%%%%%%first,second, third , and 4th nearest neighbors all in one
%%%%%%%%%%%%%plot
w=100
kl=1
figure

%subplot(1,3,1)
for j=1:4
al =[]
al1=[]
%w=w+100

kl=kl+1;
for i=1:500

 
pos = find(strcmp({pdbstruct.Model(i+1).Atom(:).AtomName},'OW'));
x =[pdbstruct.Model(i+1).Atom(pos).X];
y =[pdbstruct.Model(i+1).Atom(pos).Y];
z =[pdbstruct.Model(i+1).Atom(pos).Z];
coord= transpose([x; y; z]);

s =unique(yel((yel(:,3)<i*1000+w)&(yel(:,3)>i*1000-w),4));

s =randi([1 1019],size(s,1),size(s,2));
s=sort(s);
si = [];

for kk=1:1019;
    ap=0;
    for ct=1:size(s,1);
        if kk == s(ct)
           ap=ap+1;
        end
    end
    if ap==0;
        si = [si;kk];
    end
end

[mIdx,mD] = knnsearch(coord(s,:),coord(si,:),'K',5,'Distance',@eu);
al =[al;mD(:,kl-1)];
[mIdx,mD] = knnsearch(coord(s,:),coord(s,:),'K',5,'Distance',@eu);
al =[al;mD(:,kl)];


%save("resid_frame_"+num2str(i)+".dat",'s','-ascii')

end

hold on
[f xi]=ksdensity(al)
plot(xi,f,"linewidth",2)
%hold on
%[f xi]=ksdensity(al1)
%plot(xi,f,"LineWidth",2)
%legend(["all molecules"," jumping molecules"],"FontSize",15)
%xlabel("Min. distance to 1st jumping molecules( \AA )","FontSize",20,"Interpreter","latex")
ylabel("PDF","FontSize",20)
%title("dT ="+num2str(2*w)+" fs","fontsize",20)
end
xlabel("Min. distance to jumping molecule (\AA )","FontSize",20,"Interpreter","latex")
xlim([0 14])
ylim([0 0.6])
legend(["1st","2nd","3rd","4th"],"FontSize",15)
box on
%title("All molecules","FontSize",20)
set(groot,'defaultAxesFontSize',15)

kl=1
figure
for j=1:4
al =[]
al1=[]
%w=w+100

kl=kl+1;
for i=1:500

pos = find(strcmp({pdbstruct.Model(i+1).Atom(:).AtomName},'OW'));
x =[pdbstruct.Model(i+1).Atom(pos).X];
y =[pdbstruct.Model(i+1).Atom(pos).Y];
z =[pdbstruct.Model(i+1).Atom(pos).Z];
coord= transpose([x; y; z]);

s =unique(yel((yel(:,3)<i*1000+w)&(yel(:,3)>i*1000-w),4));

s =randi([1 1019],size(s,1),size(s,2));
s=sort(s);
[mIdx,mD] = knnsearch(coord(s,:),coord,'K',5,'Distance',@eu);
[mIdx1,mD1] = knnsearch(coord(s,:),coord(s,:),'K',5,'Distance',@eu);
al =[al;mD(:,kl)];
al1=[al1;mD1(:,kl)];

end

%hold on
%[f xi]=ksdensity(al)
%plot(xi,f,"linewidth",2,"LineStyle","--")
hold on
[f xi]=ksdensity(al1)
plot(xi,f,"LineWidth",2)
%legend(["all molecules"," jumping molecules"],"FontSize",15)
%xlabel("Min. distance to 1st jumping molecules( \AA )","FontSize",20,"Interpreter","latex")
ylabel("PDF","FontSize",20)
%title("dT ="+num2str(2*w)+" fs","fontsize",20)
end
xlabel("Min. distance to jumping molecule (\AA )","FontSize",20,"Interpreter","latex")
xlim([0.5 14])
hold on
legend(["1st","2nd","3rd","4th"],"FontSize",15)
box on

%title("Jumping Molecules","FontSize",20)
%text(0.25,0.75,"dT = "+num2str(2*w)+" fs","FontSize",20,"FontWeight","bold")

w=100
r=[]
for i=1+w:100:500000
s =size(unique(yel((yel(:,3)<i+w)&(yel(:,3)>i*1-w),4)),1);
r=[r s];
end

figure
%subplot(1,3,2)
[f xi]=ksdensity(r)
plot(xi,f,"LineWidth",3)
xlabel("Number of jumping molecules","fontsize",20)
ylabel("PDF")
%xlim([50 130])
box on

hold on
w=100
r=[]
for i=1+w:100:500000
s =size(unique(yel40((yel40(:,3)<i+w)&(yel40(:,3)>i*1-w),4)),1);
r=[r s];
end


[f xi]=ksdensity(r)
plot(xi,f,"LineWidth",3)
xlabel("Number of jumping molecules","fontsize",20)
ylabel("PDF")
%xlim([50 130])
box on



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




w=0
figure
for j=1:4
al =[]
al1=[]
w=w+100
for i=1:500

pos = find(strcmp({pdbstruct.Model(i+1).Atom(:).AtomName},'OW'));
x =[pdbstruct.Model(i+1).Atom(pos).X];
y =[pdbstruct.Model(i+1).Atom(pos).Y];
z =[pdbstruct.Model(i+1).Atom(pos).Z];
coord= transpose([x; y; z]);

s =unique(yel((yel(:,3)<i*1000+w)&(yel(:,3)>i*1000-w),4));
[mIdx,mD] = knnsearch(coord(s,:),coord,'K',2,'Distance',@eu);
[mIdx1,mD1] = knnsearch(coord(s,:),coord(s,:),'K',2,'Distance',@eu);
al =[al;mD(:,2)];
al1=[al1;mD1(:,2)];

end



subplot(2,2,j)
[f xi]=ksdensity(al)
plot(xi,f,"linewidth",2)
hold on
[f xi]=ksdensity(al1)
plot(xi,f,"LineWidth",2)
legend(["all molecules"," jumping molecules"],"FontSize",15)
xlabel("Min. distance to 1st jumping molecules( \AA )","FontSize",20,"Interpreter","latex")
ylabel("PDF")
title("dT ="+num2str(2*w)+" fs","fontsize",20)
end
%D2  = @(XI,XJ)sqrt(sum(((XI'-XJ')-31.197*round((XI'-XJ')/31.197)).^ 2)); 
%D2  = @(XI,XJ)sqrt(sum( ((XI-XJ)-31.197*round((XI-XJ)/31.197)).^ 2,2) ); 
%d1=pdist2(coord,coord,D2);
%d2=pdist2(coord,coord,@eu);



%%%2nd nearest neighbor

w=0
figure
for j=1:4
al =[]
al1=[]
w=w+100
for i=1:500

pos = find(strcmp({pdbstruct.Model(i+1).Atom(:).AtomName},'OW'));
x =[pdbstruct.Model(i+1).Atom(pos).X];
y =[pdbstruct.Model(i+1).Atom(pos).Y];
z =[pdbstruct.Model(i+1).Atom(pos).Z];
coord= transpose([x; y; z]);

s =unique(yel((yel(:,3)<i*1000+w)&(yel(:,3)>i*1000-w),4));
[mIdx,mD] = knnsearch(coord(s,:),coord,'K',5,'Distance',@eu);
[mIdx1,mD1] = knnsearch(coord(s,:),coord(s,:),'K',5,'Distance',@eu);
al =[al;mD(:,5)];
al1=[al1;mD1(:,5)];

end



subplot(2,2,j)
[f xi]=ksdensity(al)
plot(xi,f,"linewidth",2)
hold on
[f xi]=ksdensity(al1)
plot(xi,f,"LineWidth",2)
legend(["all molecules"," jumping molecules"],"FontSize",15)
xlabel("4th nearest distance to jumping molecule( \AA )","FontSize",20,"Interpreter","latex")
ylabel("PDF")
title("dT ="+num2str(2*w)+" fs","fontsize",20)
end

%%%Fourth nearest neighbor















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Jumping molecules


w=0
figure
for j=1:4
al =[]
al1=[]
w=w+100
for i=1:500

pos = find(strcmp({pdbstruct.Model(i+1).Atom(:).AtomName},'OW'));
x =[pdbstruct.Model(i+1).Atom(pos).X];
y =[pdbstruct.Model(i+1).Atom(pos).Y];
z =[pdbstruct.Model(i+1).Atom(pos).Z];
coord= transpose([x; y; z]);

s =unique(yel((yel(:,3)<i*1000+w)&(yel(:,3)>i*1000-w),4));
ss =unique(yello((yello(:,3)<i*1000+w)&(yello(:,3)>i*1000-w),4));
[mIdx,mD] = knnsearch(coord(s,:),coord,'K',3,'Distance',@eu);
[mIdx1,mD1] = knnsearch(coord(s,:),coord(s,:),'K',3,'Distance',@eu);
[mIdx2,mD2] = knnsearch(coord(ss,:),coord(ss,:),'K',3,'Distance',@eu);
al =[al;mD(:,3)];
al1=[al1;mD1(:,3)];

end

subplot(2,2,j)
[f xi]=ksdensity(al)
plot(xi,f,"linewidth",2)
hold on
[f xi]=ksdensity(al1)
plot(xi,f,"LineWidth",2)
legend(["all molecules"," jumping molecules"])
xlabel("Nearest neigbor distance( \AA )","FontSize",20,"Interpreter","latex")
ylabel("PDF")
title("dT ="+num2str(2*w),"fontsize",20)

end






