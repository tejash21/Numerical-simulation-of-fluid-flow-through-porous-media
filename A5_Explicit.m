clear
clc
t_f=input("enter value of final time= ");
x0=0;xf=1;L=xf-x0;%from lower & higher value in x find length L
y0=0;yf=1;W=yf-y0;%from lower & higher value in y find width W
Nx=input("enter value of nodes in x-direction= ");
Ny=input("enter value of nodes in y-direction= ");
dx=L/(Nx-1);%derive stepsize in x
dy=L/(Ny-1);%derive stepsize in y
a1=input("enter value of RHS of stability condition= ");%enter value less than or equal to 0.25
dt=(a1*dx*dx);%derive stepsize in time from stability condition
nt=t_f/dt;
Nt=round(nt,0);%derive max iteration in time
x=0:dx:1;%derive x axis array
y=0:dy:1;%derive y axis array
t=linspace(0,t_f,Nt);%derive time array

T=zeros(Nt,Nx,Ny);%intialize temperature array as 3D
%this is I.C.
T(1,:,:)=0;

%this is B.C.
T(2:Nt,1,:)=1;
T(2:Nt,Nx,:)=1;
T(2:Nt,:,1)=1;
T(2:Nt,:,Ny)=1; 
 
for i=2:Nt
     for j =2:Nx-1
         for k=2:Ny-1
             T(i,j,k)= T(i-1,j,k)+a1*((T(i-1,j+1,k)-2*T(i-1,j,k)+T(i-1,j-1,k))+(T(i-1,j,k+1)-2*T(i-1,j,k)+T(i-1,j,k-1)));
         end
     end
end

%now plotting
for i=1:Nx
    for j=1:Ny
        f_1(i,j)=T(Nt,i,j);
    end
end
figure(1)
pcolor(f_1);
colorbar
xlabel('Nx');
ylabel('Ny');
title('Temperature Distribution');
figure(2)
surf(f_1);
colorbar

f_3=T(:,10,10);
figure(3)
plot(t,f_3);
xlabel('time')
ylabel('Temperature for nx=10,ny=10(grid 20X20)')