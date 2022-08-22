clear
clc
L=input("enter value of length= ");
t=input("enter value of time= ");
dx=input("enter the value of given dx= ");
s=input("enter RHS number for stability condition= "); %s value should be less or equal to 0.5 
Question=input("enter the case number= "); %enter 1 for case1 & enter 2 for case2
N=L/dx; %N is nodes
dt=s*dx*dx; %This is condition to find dt.
t_steps=t/dt;
t_s=round(t_steps,0); %t_s is final step for time iteration
x=linspace(0,L,N+1); %defining nodes
t=linspace(0,t,t_s);%defining time steps

%getting values by excat solution at t=0.85
u_e = zeros(1,N+1); 
for i = 1:N+1
    u_e(1,i) = exp(t(t_s)-x(i));
end

%getting analytical value 
u=zeros(t_s,N+1);

%defining boundry condition 
for j=1:t_s
    u(j,1)=exp(t(j));
end
%defining initial condition
for k=1:N+1
    u(1,k)= exp(-x(k));
end

if Question==1
    for l=2:t_s
        for m =2:N
            u(l,m)=(u(l-1,m)*(1-dt))+((s)*(u(l-1,m+1)-2*u(l-1,m)+u(l-1,m-1)))-((0.5*dt/dx)*(u(l-1,m+1)-u(l-1,m-1)));
        end
        u(l,N+1)=u(l,N-1)*(1-2*dx); %this is boundry condition
    end
else
    for l=2:t_s
        for m =2:N
            u(l,m)=(u(l-1,m)*(1-dt))+((s)*(u(l-1,m+1)-2*u(l-1,m)+u(l-1,m-1)))-((dt/dx)*(u(l-1,m)-u(l-1,m-1)));
        end
        u(l,N+1)=u(l,N)*(1-dx); %this is boundry condition
    end
end

%plot values at final time by both exact and analytical
figure(1)
plot(x,u(end,:),'r')
hold on
plot(x,u_e,'ob')
xlabel('x')
ylabel('u(t,x)')
hold off

%now define error
e=zeros(1,N+1);
for n=1:N+1
    e(1,n)=abs(((u_e(1,n)-u(end,n))*100)/u(end,n));
end

%plot error
figure(2)
plot(x,e)
xlabel('x')
ylabel('e(%)')