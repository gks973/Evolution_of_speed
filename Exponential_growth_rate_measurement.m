clc; clear all;

dt=0.001;   %timestep
N=1;        %population size


beta=1;
c=1;  %ratio between latency time and infectious period duration (c = tau/T)
k=10; %shape factor of the Gamma distribution  of latency times and infectious periods

if c==0 %c cannot be exactly 0, as this results in division by 0
    c=0.01;
end

time_max=10; %what is the maximal infectious period duration that we simulate?
dtime=0.05;  %by which stepsize should we vary infectious period duration?


growth_rate=zeros(1,10);

for i=1:time_max/dtime %how many infectious period durations should we simulate?
    
    time=0.05*i; % infectious period duration T
    tau=time*c;  % exposed state duration / latency time
    gamma=1/time;% recovery rate
    eta=1/tau;   % rate of progression from exposed to infectious state
    T=ceil(2*(1/gamma+1/eta)); % total simulation time is equal to two full disease generations
    
    if tau/k<10*dt||time/k<10*dt %if a subcompartment becomes shorter than 10*dt, decrease dt by a factor 10
        dt=dt/10;
    else
        dt=0.001;
    end
    
    S=zeros(1,round(T/dt));
    E=zeros(k,round(T/dt));
    I=zeros(k,round(T/dt));
    R=zeros(1,round(T/dt));
    
    growth_over_time=zeros(1,round(T/dt));

    S(1)=1;
    E(:,1)=10^-12/k;
    I(:,1)=10^-12/k; %start out the epidemic with 10^-12 exposed and 10^-12 infectious
                     %equally distributed in all k subcompartments
    dE=zeros(k,1);
    dI=zeros(k,1);
    
    for t=1:T/dt
        
        dS=-beta*S(t)*sum(I(:,t))/N;                   %solve differential equations
        dE(1)=beta*S(t)*sum(I(:,t))/N-(k*eta).*E(1,t); %there are k E- and I-compartments
        for n=2:k                                      %in order to obtain a Gamma-distributed
            dE(n)=(k*eta)*E(n-1,t)-(k*eta).*E(n,t);    %duration with shape factor k
        end
        dI(1)=k*eta*E(k,t)-(k*gamma).*I(1,t);
        for n=2:k
            dI(n)=(k*gamma).*I(n-1,t)-(k*gamma).*I(n,t);
        end
        dR=k*gamma.*I(k,t);

        S(t+1)=S(t)+dS*dt;
        E(:,t+1)=E(:,t)+dE(:)*dt;
        I(:,t+1)=I(:,t)+dI(:)*dt;
        R(t+1)=R(t)+dR*dt;
        
        if round(1*(time+tau)/dt)<t&&t<round(2*(time+tau)/dt)     %calculate the epidemic growth rate over the
            growth_over_time(t)=(sum(I(:,t))/sum(I(:,t-1))-1)/dt; %second disease generation
        end
        
        if sum(I(:,t))<0 %stop the computation if I becomes negative due to numerical error
            break
        end
    end
    
    
    growth_rate(i)=mean(growth_over_time(round(1*(time+tau)/dt):round(2*(time+tau)/dt)));
    %Calculate the mean of the growth rate to smooth out fluctuations.
    %Fluctuations can be large over the course of a disease generation when
    %k and c are large
    
end

[max_rate,max_point]=max(growth_rate); %find location and value of the maxima of the growth rate

plot(dtime:dtime:time_max,growth_rate)
xlabel('c')
ylabel('Exponential growth rate')
ylim([0 0.3])
