%________________________________________________________________________%
% Diversity-enhanced fuzzy Multi-Objective Particle Swarm Optimizer      %
% (f-MOPSO/Div)                                                          %
%                                                                        %
% This code is usable to solve multi/many-objective optimization problems%
%                                                                        %
% Developed in MATLAB R2018b                                             %
%                                                                        %
% Author and programmer: Farshad Rezaei, PhD                             %
%                                                                        %
% e-Mail: farshad.rezaei@gmail.com                                       %
%         f.rezaei@alumni.iut.ac.ir                                      %
%                                                                        %
% Homepage: https://www.linkedin.com/in/farshad-rezaei-5a92559a/         %
%                                                                        %
% Main paper: Rezaei, F., Safavi, H.R.,(2020) "f-MOPSO/Div: an improved  %
% extreme-point-based multi-objective PSO algorithm applied to a         %
% socio-economic-environmental conjunctive water use problem",           %
% Environmental Monitoring and Assessment, 192(12): 1-27                 %
%________________________________________________________________________%
% f-MOPSO/Div algorithm 
clc
clear
close all
tic
run=1; % Maximum number of the algorithm runnings conducted
invertgen=zeros(run);
maxit=1000; % Maximum number of the iterations
stallit=1000; % Stall iterations: Number of the iterations over which the global guides are not getting better
np=50; % Population size
Function_name='DTLZ6'; % Name of the test function selected here from the DTLZ problems family
[lb,ub,nx,nobj,fobj]=Objective_Function(Function_name); % Load details of the selected benchmark function
varmin=lb; % Upper bound defined for the positions which can generally be a desired vector
varmax=ub; % Lower bound defined for the positions which can generally be a desired vector
limvel=0.1; % A ratio of the maximum distance in the search space resulting in the maximum velocity
velmax=limvel*(varmax(1,1:nx)-varmin(1,1:nx)); % Upper bound defined for the velocities
velmin=-velmax; % Lower bound defined for the velocities
st=0.5; % Weight step as a constant parameter in this algorithm
stallpower=6; % Power of 10 denoting the maximum allowable error between the objective values over the stall iterations
indexx_main=zeros(run,maxit*np); % An index calculated to delineate if a solution is non-dominated
nond_archive_main=zeros(run,maxit*np,nx); % Elements of a solution stored in the non-dominated archive
znp_main=zeros(run,nobj,maxit*np); % The objectives calculated for a solution stored in the non-dominated archive
fitness_main=zeros(run,maxit); % Average Dominance Indices of the population at each iteration
obj_min=zeros(nobj); % Minimum value found for each objective among the archived non-dominated solutions
t_best=zeros(nobj); % Indicator of minimum value of each objective in the non-dominated archive
pareto_size=1000; % Size of the true Pareto-set
obj=zeros(nobj,max(maxit*np,pareto_size)); % Objectives of the non-dominated solutions stored in the archive
ideal=zeros(nobj,1); % Coordinates of the ideal point of the final Pareto-front
distance=zeros(pareto_size); % Pre-allocating distance for speed
a=zeros(max(nobj,pareto_size),max(nobj,pareto_size)); % Pre-allocating a for speed
pareto=zeros(maxit*np,nx); % Elements (decision variables) of the non-dominated archived solutions
d=zeros(1,maxit*np); % Pre-allocating d for speed
k_max=0.9; % maximum value of k used to tune the constriction coefficient at each iteration
k_min=0.4; % minimum value of k used to tune the constriction coefficient at each iteration
mean_ideal_dis=zeros(run); % Pre-allocating for speed
spacing=zeros(run); % Pre-allocating for speed
no_solutions=zeros(run); % Pre-allocating for speed
ideal_point=zeros(run,nobj); % Pre-allocating for speed
for nrun=1:run
    [indexx,nond_archive,znp,fitness,np_archive,maxit_final]=f_MOPSO_Div(maxit,np,nx,nobj,st,k_max,k_min,varmin,varmax,velmin,velmax,stallit,stallpower,fobj);
    indexx_main(nrun,1:maxit*np)=indexx(1,1:maxit*np);
    nond_archive_main(nrun,:,1:nx)=nond_archive(:,1:nx);
    znp_main(nrun,:,:)=znp(:,:);
    fitness_main(nrun,:)=fitness(1,:);
    x1=zeros(maxit_final);
    y1=zeros(maxit_final);
    d1=0;
    for yy=1:np_archive
        if indexx_main(nrun,yy)==1
            d1=d1+1;
            pareto(d1,:)=nond_archive_main(nrun,yy,:);
            for i=1:nobj
                obj(i,d1)=znp_main(nrun,i,yy);
            end
            if d1==1
                for i=1:nobj
                    obj_min(i)=obj(i,d1);
                    t_best(i)=d1;
                end
            else
                for i=1:nobj
                    if obj(i,d1)<obj_min(i)
                        obj_min(i)=obj(i,d1);
                        t_best(i)=d1;
                    end
                end
            end
        end
    end
    if nobj==2
        obj1=obj(1,1:d1);
        obj2=obj(2,1:d1);
        figure(nrun);
        plot(obj1,obj2,'ro')
        xlabel('Z1');
        ylabel('Z2');
        hold on
    elseif nobj==3
        obj1=obj(1,1:d1);
        obj2=obj(2,1:d1);
        obj3=obj(3,1:d1);
        figure(nrun);
        plot3(obj1,obj2,obj3,'ro')
        xlabel('Z1');
        ylabel('Z2');
        zlabel('Z3');
        hold on
    end
    for i=1:nobj
        ideal(i,1)=min(obj(i,1:d1));
    end
    
%     disp('The Pareto = ');
%     for i = 1:d1
%         for j = 1:nx
%             disp(num2str(pareto(i,j)));
%         end
%     end

%     id=ideal(:,1); % For a problem whose true ideal point is unknown
    id=zeros(nobj,1); % For a problem whose true ideal point is known such as that in the DTLZ family
    sum2=0;
    for i=1:d1
        sum1=0;
        for j=1:nobj
%             sum1=sum1+((obj(j,i))-ideal(j,1))^2;
            sum1=sum1+((obj(j,i))-id(j,1))^2;
        end
        sum2=sum2+sqrt(sum1);
    end
    mid=sum2/d1;
    for i=1:d1
        min_ss=inf;
        for j=1:d1
            if i~=j
                ss=0;
                for m=1:nobj
                    ss=ss+abs(obj(m,i)-obj(m,j));
                end
                if ss<min_ss
                    min_ss=ss;
                end
            end
        end
        d(1,i)=min_ss;
    end
    d_average=mean(d(1,1:d1));
    sum2=0;
    for i=1:d1
        sum2=sum2+(d_average-d(1,i))^2;
    end
    s=sqrt(sum2/(d1-1));
    ns=d1;

    mean_ideal_dis(nrun)=mid;
    spacing(nrun)=s;
    no_solutions(nrun)=ns;
    ideal_point(nrun,:)=ideal(:);
end
disp('Performance Metrics in total runs =');

% Computing the Performance Criteria

% MID = Mean Ideal Distance
disp(['Average MID = ',num2str(mean(mean_ideal_dis(1:run)))]); 
disp(['Best MID = ',num2str(min(mean_ideal_dis(1:run)))]); 
disp(['Worst MID = ',num2str(max(mean_ideal_dis(1:run)))]); 
disp(['Std of MID = ',num2str(std(mean_ideal_dis(1:run)))]); 

% S = Spacing
disp(['Average S = ',num2str(mean(spacing(1:run)))]); 
disp(['Best S = ',num2str(min(spacing(1:run)))]); 
disp(['Worst S = ',num2str(max(spacing(1:run)))]); 
disp(['Std of S = ',num2str(std(spacing(1:run)))]); 

% NS = Number of Solutions
disp(['Average NS = ',num2str(mean(no_solutions(1:run)))]); 
disp(['Best NS = ',num2str(max(no_solutions(1:run)))]); 
disp(['Worst NS = ',num2str(min(no_solutions(1:run)))]); 
disp(['Std of NS = ',num2str(std(no_solutions(1:run)))]); 

% Ideal Point
for i=1:nobj
    disp(['Average Minimum ',num2str(i),'th Objective = ',num2str(mean(ideal_point(1:run,i)))]);
end
toc