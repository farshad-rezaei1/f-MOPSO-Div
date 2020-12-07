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
function [indexx,nond_archive,znp,fitness,np_archive,maxit_final]=f_MOPSO_Div(maxit,np,nx,nobj,st,k_max,k_min,varmin,varmax,velmin,velmax,stallit,stallpower,fobj) 

% MOPSO Parameters
maxit_final=maxit;
stall=zeros(nobj);
pp_pbest=zeros(np,nx);
nond_archive=zeros(maxit*np,nx);
dom_archive=zeros(maxit);
fitness=zeros(maxit);
z_archive=zeros(nobj,maxit);
z_ppbest=zeros(nobj,np);
znp=zeros(nobj,maxit*np);
zobj=zeros(nobj);
z_min=inf*ones(nobj);
z_max=-inf*ones(nobj);
p_max=1; % Maximum dynamic mutation probability 
p_min=1/nx; % Minimum dynamic mutation probability
nm=20;
indexx=ones(1,maxit*np);

it=1;
np_archive=0;

% Initialization process of the algorithm
[pp,pv]=Initialization(np,nx,varmax,varmin,velmax,velmin);
for j=1:np
    np_archive=np_archive+1;
    x(1:nx)=pp(j,1:nx); 
    
    % Function Evaluation   
    [zobj]=fobj(x);
    z_ppbest(1:nobj,j)=zobj(1:nobj);
    for m = 1:nobj
        if zobj(m)<z_min(m)
            z_min(m)=zobj(m);
        elseif zobj(m)>z_max(m)
            z_max(m)=zobj(m);
        end
    end
    for i=1:nobj
        znp(i,np_archive)=zobj(i);
    end
    nond_archive(np_archive,1:nx)=x(1:nx);
    if np_archive~=1
        [indexx]=Archive(indexx,nobj,znp,np_archive);
    end
end

% SFIS Implementation
[dom,dom_ppbest]=SFIS(np_archive,nobj,zobj,z_ppbest,znp,z_min,z_max,st,it);
for j=1:np
    dom_ppbest(j)=dom(j);
    pp_pbest(j,1:nx)=pp(j,1:nx);
    for m=1:nobj
        z_ppbest(m,j)=zobj(m);
    end
end

% Determining the Global Best Particle
[g_elem]=Gbest_Selection(z_ppbest);
pp_gbest(1:nx)=pp_pbest(g_elem,1:nx);
dom_gbest=dom_ppbest(g_elem);
 
% Sending the Gbest Particle to an Archive
for i=1:nobj
    z_archive(i,it)=z_ppbest(i,g_elem);
end
dom_archive(it)=dom_gbest;
sumf=0;
for j=1:np
    sumf=sumf+dom(j);
end
fitness(1,it)=sumf/np;

% Main Loop

while it<maxit
    it = it+1;
    k=k_max-(k_max-k_min)*(it-2)/(maxit-2); 
%     disp(['Number of Iterations= ',num2str(it)]);
    for j = 1:np
        phi1=2.05*rand(1,nx);
        phi2=2.05*rand(1,nx);
        phi=phi1+phi2;
        khi=(2*k)./abs(2-phi-sqrt(phi.*(phi-4))); 
        
        % Updating the velocity of the particles
        pv(j,1:nx)=khi.*(pv(j,1:nx)+phi1.*(pp_pbest(j,1:nx)-pp(j,1:nx))+phi2.*(pp_gbest(1:nx)-pp(j,1:nx)));
        
        % Returning the velocity of the particles if going beyond the velocity boundaries
        flag4lbv=pv(j,:)<velmin(1,:);
        flag4ubv=pv(j,:)>velmax(1,:);
        pv(j,:)=(pv(j,:)).*(~(flag4lbv+flag4ubv))+velmin.*flag4lbv+velmax.*flag4ubv;
        
        % Updating the position of the particles
        pp(j,:)=pp(j,:)+pv(j,:);
        
        % Imposing Mutation
        for i=1:nx
            rm=rand(1,1);
            pm=p_max-((p_max-p_min)/maxit)*it;
            if rm<pm
                r=rand(1,1);
                if r<0.5
                    delta=(2*r)^(1/(nm + 1))-1;
                    pp(j,i)=pp(j,i)+delta*(varmax(1,i)-varmin(1,i));
                else
                    delta=1-((2*(1 - r))^(1/(nm + 1)));
                    pp(j,i)=pp(j,i)+delta*(varmax(1,i)-varmin(1,i));
                end    
            end
        end
        
        % Returning the position and velocity of the particles if going beyond the position boundaries
        flag4lbp=pp(j,:)<varmin(1,:);
        flag4ubp=pp(j,:)>varmax(1,:);
        pp(j,:)=(pp(j,:)).*(~(flag4lbp+flag4ubp))+varmin.*flag4lbp+varmax.*flag4ubp; 
        pv(j,:)=(pv(j,:)).*(ones(1,nx)-2*(flag4lbp+flag4ubp));        
        np_archive=np_archive+1;
        x(1:nx)=pp(j,1:nx);
        
        % Function Evaluation   
        [zobj]=fobj(x);
        for m=1:nobj
            if zobj(m)<z_min(m)
                z_min(m)=zobj(m);
            elseif zobj(m)>z_max(m)
                z_max(m)=zobj(m);
            end
        end
        for i=1:nobj
            znp(i,np_archive) = zobj(i);
        end
        nond_archive(np_archive,1:nx)=x(1:nx);
        if np_archive~=1
            [indexx]=Archive(indexx,nobj,znp,np_archive);
        end
    end
    
    % SFIS Implementation
    [dom,dom_ppbest]=SFIS(np_archive,nobj,zobj,z_ppbest,znp,z_min,z_max,st,it);
    for j=1:np
        if dom(j)<dom_ppbest(j)
            dom_ppbest(j)=dom(j);
            pp_pbest(j,1:nx)=pp(j,1:nx);
            for m=1:nobj
                z_ppbest(m,j)=zobj(m);
            end
        end
    end
    
    % Determining the Global Best Particle
    [g_elem]=Gbest_Selection(z_ppbest);
    pp_gbest(1:nx)=pp_pbest(g_elem,1:nx);
    dom_gbest=dom_ppbest(g_elem);

    % Sending the Gbest Particle to an Archive
    for i=1:nobj
        z_archive(i,it)=z_ppbest(i,g_elem);
    end
    dom_archive(it)=dom_gbest;
    sumf=0;
    for j=1:np
        sumf=sumf+dom(j);
    end
    fitness(1,it)=sumf/np;
   
    % Termination Criterion
    if it>=stallit
        stall_no=0;
        for i=1:nobj
            stall(i)=z_archive(i,it)-z_archive(i,it-stallit+1);
            if abs(stall(i))<10^(-stallpower)
                stall_no=stall_no + 1;
            end
        end
        if stall_no==nobj
            maxit_final=it;
            break
        end
    end
end