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
function [g_elem]=Gbest_Selection(z_ppbest)
nobj=size(z_ppbest,1);
np=size(z_ppbest,2);
part=zeros(np);
z_div=zeros(nobj,np);
flag=ones(1,np);
div=zeros(np);

% Preserving the Pbest Particles Once Detected to be Non-dominated
for j=1:np
    if j~=1
        for y=1:j-1
            if flag(1,y)==1
                c1=0;c2=0;
                for i=1:nobj
                    if (z_ppbest(i,j)<=z_ppbest(i,y)) 
                        c1=c1+1;
                    elseif (z_ppbest(i,j)>=z_ppbest(i,y))
                        c2=c2+1;
                    end
                end
                if c1==nobj 
                    flag(1,y)=0;
                elseif c2 == nobj && c2 ~= c1
                    flag(1,j)=0;
                    break
                end
            end
        end
    end
end

% Descriminating the Non-dominated Pbests from the whole Pbests' Population
d=0;
for j=1:np
    if flag(1,j)==1
        d=d+1;
        part(d)=j;
        for i=1:nobj
            z_div(i,d)=z_ppbest(i,j);
        end
    end
end

% Evaluating the Diversity of the Pbests
for i=1:d
    min_ss=inf;
    for j=1:d
        ss=0;
        for m=1:nobj
            ss=ss+abs(z_div(m,i)-z_div(m,j));
        end
        if i~=j && ss<min_ss
            min_ss=ss;
        end
    end
    div(i)=min_ss;
end

% Selecting the Most-diversified Pbest as the Gbest Particle
max_div=-inf;
for i=1:d
    if div(i)>=max_div
        max_div=div(i);
        t=part(i);
    end
end
g_elem=t; % g_elem is the position of the Gbest particle in the current population 
end