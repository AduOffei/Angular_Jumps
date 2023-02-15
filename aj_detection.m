% In this program I extract angular magnitudes and times for the molecules
% here import the dipole/hh vector
A = importdata("dp_space.dat");
t_step=1
num_mol=1019
fc = 1/100;

tim = size(A,1)

% the loop over particles begins here
b_time = [];
b_ang  = [];
b_dt   = [];
b_mol  = [];
for j  = 1:num_mol
%I do butterworth filtering here to obtain fxyz a filerd vector
xyz=A(:,1+3*(j-1):3+3*(j-1));

fs = 1/t_step;
[b,a] = butter(2,fc/(fs/2));
%freqz(b,a);
fxyz = filtfilt(b,a,xyz);
sxyz = zeros(size(fxyz));
% mean filtering here
sxyz(:,1) = smooth(fxyz(:,1),100/t_step);
sxyz(:,2) = smooth(fxyz(:,2),100/t_step);
sxyz(:,3) = smooth(fxyz(:,3),100/t_step);

for i=1:tim
    sxyz(i,:)=sxyz(i,:)/sqrt(dot(sxyz(i,:),sxyz(i,:)));
end

% computing the velocity and dot product
vxyz =diff(sxyz);
qxyz = cross(sxyz(1:end-1,:),vxyz);


for i=1:tim-1
    qxyz(i,:)=qxyz(i,:)/sqrt(dot(qxyz(i,:),qxyz(i,:)));
end
% computing the dot product between consecutive frames as indication of
% change of geodesic

ah = zeros(size(qxyz,1)-1,1);
for i= 1:size(qxyz,1)-1
    ah(i)=1.0000- dot(qxyz(i,:),qxyz(i+1,:))/sqrt(dot(qxyz(i,:),qxyz(i,:))*dot(qxyz(i+1,:),qxyz(i+1,:)));
    if ah(i) <=10^-4
        ah(i)=0.0;
    end
end

l = islocalmax(ah);

% here I generate an array of l indicating the presence of a jump t with
t = find(l ==1);

% change in jump of angle , change in the angle
 ang = zeros(size(t));
 dt = zeros(size(t));
 mol = zeros(size(t));
 ang(1) = dot(  xyz(t(1),:), xyz(1,:))/sqrt(dot(xyz(t(1),:),xyz(t(1),:))*dot(xyz(1,:),xyz(1,:)));
 ang(1) = acos(ang(1));
 dt(1) = t(1); 
 mol(:) = j;
 
for i = 2:size(t,1)
 ang(i) = dot( xyz(t(i),:),xyz(t(i-1),:) )/sqrt( dot( xyz(t(i),:),xyz(t(i),:) ) * dot(xyz(t(i-1),:),xyz(t(i-1),:))  );
 ang(i) = acos(ang(i));
 dt(i) = t(i)-t(i-1);
end
ang= ang*(180/pi);
b_ang = [b_ang;ang];
b_dt =[b_dt;dt];
b_time =[b_time;t];
b_mol = [b_mol;mol];
end





tbs =[b_ang,t_step*b_dt,t_step*b_time,b_mol];
[~,idx] = sort(tbs(:,3)); % sort just the first column
cbs = tbs(idx,:);   % sort the whole matrix using the sort indices
save("gd_angle_vs_time_dp_1.mat",'cbs','-v7.3');
save("gd_angle_vs_time_dp_1.dat",'cbs','-ascii');






