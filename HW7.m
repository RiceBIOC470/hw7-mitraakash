% Akash Mitra
% am132

%HW7

% Problem 1: Modeling population growth
% The simplest model for a growing population assumes that each current
% individual has equal likelihood to divide, which yields a differential
% equation dx/dt = a*x where a is the division rate. It is easy to see that
% this will grow exponentially without bound. A simple modification is to
% assume that the growth rate slows done as the population reaches some
% maximum value N so that dx/dt = a*x*(1-x/N). Defining X = x/N, we have 
% dX/dt = a*X*(1-X).  
% Part 1. This equation has two fixed points at 0 and 1. Explain the
% meaning of these two points.

% If the fixed points refer to the dx/dt variable the,
% Two fixed points of 0 and 1 exist to denote the sigmoidal shape of the
% curve. At a value of 0, the delta change of growth in population is 0 and
% the curve is saturated. The growth rate slows, and eventually stops. The
% delta value of 1 implies that growth is increasing exponentially, when 
% the population/resources is not saturated. 

% If the fixed points refer to the value of x, the fixed points of 0 and 1
% would have no population growth. This is due to the fact that at x=0,
% there is no population to grow, and at x=1, it would also not be possible
% to increase the population.

% Part 2: Evaluate the stability of these fixed points. Does it depend on
% the value of the parameter a? 

% In this situation, x represents the population size, N the carrying
% capacity and a defines the growth rate. Early, unimpeded growth rate is
% modeled by the first term +aX. The value of the rate a represents the
% proportional increase of the population x in one unit of time. Later, as
% the population grows, the modulus of the second term becomes almost as 
% large as the first, as some members of the population x interfere with 
% each other by competing for some critical resource

% In this situation, the a defines the growth rate but since there is no
% population to expand, the population remains stable at x=0 and x=1

% Part 3: Write a function that takes two inputs - the initial condition x0
% and the a parameter and integrates the equation forward in time. Make
% your code return two variables - the timecourse of X and the time
% required for the population to reach 99% of its maximum value. 


t=[];
j=1;
    while (t(j) < 0.99)
        rhs = @(t,x,a,N) a*x*(1-x/N);
        dt = 0.1;
        interval = [0 100];
        nstep = (interval(2) - interval(1))/dt;
        sol1(1) = 0;
        for i = 2:nstep
            sol1(i) = sol1(i-1) + rhs(i,10,100,1000) * dt;
        end
        t = linspace(interval(1), interval(2), nstep);
        j=j+1;
    end
disp(j);





% Part 4: Another possible model is to consider discrete generations
% instead allowing the population to vary continuously. e.g. X(t+1) = a*
% X(t)*(1-X(t)). Consider this model and vary the a parameter in the range 0
% < a <= 4. For each value of a choose 200 random starting points  0 < x0 < 1 
% and iterate the equation forward to steady state. For each final
% value Xf, plot the point in the plane (a,Xf) so that at the end you will
% have produced a bifucation diagram showing all possible final values of
% Xf at each value of a. Explain your results.
        
  
for a=0.1:0.1:4
    for i=1:200
        x = rand(1,1);
        xf = x;
        for i = 1:10
            xf = a*xf*(1-xf);
        end
        plot(a,xf,'.');
        hold on;
    end
end

% The bifurcation plot may arise due to the fact that the population growth
% is being regulated in different states. There is an initial low steady
% state, followed by an area of instability, and eventually moves towards a
% high state of stability.

% Problem 2. Genetic toggle switches. 
% Consider a genetic system of two genes A and B in which each gene
% product represses the expression of the other. Make the following
% assumptions: 
% a. Repression is cooperative:  each promotor region of one gene has 4
% binding sites for the other protein and that all of these need to be
% occupied for the gene to be repressed. 
% b. You can ignore the intermediate mRNA production so that the product of
% the synthesis of one gene can be assumed to directly repress the other
% c. the system is prefectly symmetric so that the degradation
% times, binding strengths etc are the same for both genes. 
% d. You can choose time and concentration scales so that all Michaelis
% binding constants and degradation times are equal to 1. 
%
% Part 1. Write down a two equation model (one for each gene product) for
% this system. Your model should have one free parameter corresponding to the
% maximum rate of expression of the gene, call it V. 

% Test condition of repressor for second gene starting high and values of 
% 1 for gene A and 2 for gene B.

% Model created in toggle_switch function

[t2,P]=ode45(@toggle_switch, 0:1:100, [3, 2, 0, 0]);


figure()
plot(t2,P(:,4),'*b','LineWidth',2); 
hold on;
plot(t2,P(:,2),'*r','LineWidth',2);
xlabel('Time'); 
ylabel('mRNA Expression');
legend('Gene A','Gene B');
set(gca,'FontSize', 18);

%
% Part 2. Write code to integrate your model in time and plot the results for V = 5 for two cases, 
% one in which A0 > B0 and one in which B0 > A0. 

[t3,P3]=ode45(@toggle_switch2, 0:1:100, [10, 3, 1, 1.5]);

figure()
plot(t3,P3(:,4),'*b','LineWidth',2); 
hold on;
plot(t3,P3(:,2),'*r','LineWidth',2);
xlabel('Time'); 
ylabel('mRNA Expression');
legend('Gene A','Gene B');
set(gca,'FontSize', 18);

[t4,P4]=ode45(@toggle_switch2, 0:1:100, [1, 3, 10, 1.5]);
figure()
plot(t4,P4(:,4),'*b','LineWidth',2); 
hold on;
plot(t4,P4(:,2),'*r','LineWidth',2);
xlabel('Time'); 
ylabel('mRNA Expression');
legend('Gene A','Gene B');
set(gca,'FontSize', 18);

%
% Part 3. By any means you want, write code to produce a bifurcation diagram showing all
% fixed points of the system as a function of the V parameter. 

randV = (100-1).*rand(100,1)+1;
figure();
for i = 1:1:20
    [t4,P4]=ode45(@toggle_switch3, 0:1:100, [5, 3, 7, 1.5, randV(i)]);
    hold on;
    plot(randV(i),P4(:,4),'*b','LineWidth',2); 
    hold on;
    plot(randV(i),P4(:,2),'*r','LineWidth',2);
    xlabel('V'); 
    ylabel('mRNA Expression');
    legend('Gene A','Gene B');
    set(gca,'FontSize', 18);
end

figure;
hold on;
ku=0;
for kb = 0:1:10
    polycoeff = [1 -kb 1 -ku];
    [t4,P4]=ode45(@toggle_switch3, 0:1:100, [5, 3, 7, 1.5, randV(i)]);
    rts
    P3(:,4) = P3(:,4)(imag(P3(:,4)) == 0);
    P3(:,2) = P3(:,2)(imag(P3(:,2)) == 0);
    plot(kb*ones(length(P3(:,4)),1,P3(:,4), 'r.'));
    plot(kb*ones(length(P3(:,2)),1,P3(:,2), 'r.'))
end
hold off;

ku = 0;
Repressor1 = 10;
Kd = 10;
n1= 4;
figure;
hold on;
func2 = @(x,V) (V*x^n1)/(1+(x^n1))-x;
for V = 0:0.05:5
    gene2 = @(x) func2(x,V);
    for in = 1:0.1:5
        [rt,~,exitflag] = fzero(gene2,in);
        if exitflag == 1
            plot(V,rt,'k.');
            hold on;
        end
    end
end
hold off;


