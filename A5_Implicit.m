clear
clc
t_f=input("enter value of final time= ");
x0=0;xf=1;L=xf-x0;%from lower & higher value in x find length L
y0=0;yf=1;W=yf-y0;%from lower & higher value in y find width W
Nx=input("enter value of nodes in x-direction= ");
Ny=input("enter value of nodes in y-direction= ");
dx=L/(Nx-1);%derive stepsize in x
dy=L/(Ny-1);%derive stepsize in y
a1=input("enter value of RHS of stability condition= ");%enter value same as in explicit to compare
dt=(a1*dx*dx);%derive stepsize in time from stability condition
nt=t_f/dt;
Nt=round(nt,0);%derive max iteration in time
x=0:dx:1;%derive x axis array
y=0:dy:1;%derive y axis array
t=linspace(0,t_f,Nt);%derive time array

%derive constants 
a1=dt/(dx*dx);
a2=dt/(dy*dy);
a3=((-2*a1)-(2*a2)-1);

%make coefficient matrix
c=zeros((Nx-2)*(Ny-2),(Nx-2)*(Ny-2));
for i=1:(Nx-2)*(Ny-2)
    c(i,i)=a3;
end
for i=1:(Nx-2)*(Ny-2)-1
    c(i,i+1)=a2;
end
for i =(Nx-2):(Nx-2):(Nx-2)*(Ny-2)-1
    c(i,i+1)=0;
end
for i=2:(Nx-2)*(Ny-2)
     c(i,i-1)=a2;
end
for i =(Nx-1):(Nx-2):(Nx-2)*(Ny-2)
    c(i,i-1)=0;
end
for i=1:(Nx-2)*(Ny-2)-(Nx-2)
    c(i,i+(Nx-2))=a1;
end
for i=(Nx-1):(Nx-2)*(Ny-2)
    c(i,i-(Nx-2))=a1;
end     

%make RHS matrix
B=zeros((Nx-2),(Ny-2));
B(1,:)=-a1;
B(Nx-2,:)=-a1;
B(:,1)=-a1;
B(:,Nx-2)=-a1;
B(1,1)=-(a1+a2);
B(1,Ny-2)=-(a1+a2);
B((Nx-2),1)=-(a1+a2);
B((Nx-2),(Ny-2))=-(a1+a2);
b=reshape(B,[],1);

f3 = zeros(Nt,1);
T=zeros((Nx-2)*(Ny-2),1);
V=T+b;

for i=1:Nt
    T_f=c\V;
    V=T_f+b; 
    Tf = reshape(T_f,(Nx-2),(Ny-2));
    f3(i,1)=Tf(9,9);
end

% now plotting
s=reshape(T_f,(Nx-2),(Ny-2));
f = zeros(Nx,Ny);
f(1,:) = 1;
f(Nx,:) = 1;
f(:,1) = 1;
f(:,Ny) = 1;
for i=2:Nx-1
    for j = 2:Ny-1
        f(i,j) = s(i-1,j-1);
    end
end

figure(1)
pcolor(f);
colorbar
xlabel('Nx');
ylabel('Ny');
title('Temperature Distribution');
figure(2)
surf(f);
colorbar

figure(3)
plot(t,f3);
xlabel('time')
ylabel('Temperature for nx=10,ny=10(grid 20X20)')


