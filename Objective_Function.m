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
function [lb,ub,nx,nobj,fobj]=Objective_Function(F)
% Insert the characteristics of your desired problem suites for function
% evaluation

switch F
    case 'DTLZ1'
        fobj=@DTLZ1;
        nobj=5;
        kk=5;
        nx=nobj+kk-1;
        lb=0*ones(1,nx);
        ub=1*ones(1,nx);
        
    case 'DTLZ2'
        fobj=@DTLZ2;
        nobj=5;
        kk=10;
        nx=nobj+kk-1;
        lb=0*ones(1,nx);
        ub=1*ones(1,nx);
        
    case 'DTLZ3'
        fobj=@DTLZ3;
        nobj=5;
        kk=10;
        nx=nobj+kk-1;
        lb=0*ones(1,nx);
        ub=1*ones(1,nx);
        
    case 'DTLZ4'
        fobj=@DTLZ4;
        nobj=5;
        kk=10;
        nx=nobj+kk-1;
        lb=0*ones(1,nx);
        ub=1*ones(1,nx);
              
    case 'DTLZ5'
        fobj=@DTLZ5;
        nobj=5;
        kk=10;
        nx=nobj+kk-1;
        lb=0*ones(1,nx);
        ub=1*ones(1,nx);
        
    case 'DTLZ6'
        fobj=@DTLZ6;
        nobj=3;
        kk=10;
        nx=nobj+kk-1;
        lb=0*ones(1,nx);
        ub=1*ones(1,nx);
end

% DTLZ1

function [zobj]=DTLZ1(x)
sum3=0;
    for i=nobj-1+1:nobj-1+kk
        sum3=sum3+(x(i)-0.5)^2-cos(20*pi*(x(i)-0.5));
    end
    m=0;
    for n=nobj:-1:1
        m=m + 1;
        zobj(m)=1;
        if n~=1
            for i=1:n-1
                zobj(m)=zobj(m)*x(i);
            end
            if n~=nobj
               zobj(m)=zobj(m)*(1-x(n));
            else
                continue
            end
        else
            zobj(m)=zobj(m)*(1-x(n));
        end
        zobj(m)=zobj(m)*(1+100*(kk + sum3))*1/2;
    end
end

% DTLZ2

function [zobj] = DTLZ2(x)
sum3=0;
    for i=nobj-1+1:nobj-1+kk
        sum3=sum3+(x(i)-0.5)^2;
    end
    m=0;
    for n=nobj:-1:1
        m=m+1;
        zobj(m)=1;
        if n~=1
            for i=1:n-1
                zobj(m)=zobj(m)*cos(x(i)*pi/2);
            end
            if n~=nobj
               zobj(m)=zobj(m)*sin(x(i)*pi/2);
            else
                continue
            end
        else
            zobj(m)=zobj(m)*sin(x(n)*pi/2);
        end
        zobj(m)=zobj(m)*(1+sum3);
    end
end

% DTLZ3

function [zobj]=DTLZ3(x)
sum3=0;
    for i=nobj-1+1:nobj-1+kk
        sum3=sum3+(x(i)-0.5)^2-cos(20*pi*(x(i)-0.5));
    end
    m=0;
    for n=nobj:-1:1
        m=m + 1;
        zobj(m)=1;
        if n~=1
            for i=1:n-1
                zobj(m)=zobj(m)*cos(x(i)*pi/2);
            end
            if n~=nobj
               zobj(m)=zobj(m)*sin(x(n)*pi/2);
            else
                continue
            end
        else
            zobj(m)=zobj(m)*sin(x(n)*pi/2);
        end
        zobj(m)=zobj(m)*(1+100*(kk + sum3));
    end
end

% DTLZ4

function [zobj]=DTLZ4(x)
alpha=100;
sum3=0;
    for i=nobj-1+1:nobj-1+kk
        sum3=sum3+(x(i)-0.5)^2;
    end
    m=0;
    for n=nobj:-1:1
        m=m+1;
        zobj(m)=1;
        if n~=1
            for i=1:n-1
                zobj(m)=zobj(m)*cos(x(i)^alpha*pi/2);
            end
            if n~=nobj
               zobj(m)=zobj(m)*sin(x(n)^alpha*pi/2);
            else
                continue
            end
        else
            zobj(m)=zobj(m)*sin(x(n)^alpha*pi/2);
        end
        zobj(m)=zobj(m)*(1+sum3);
    end
end

% DTLZ5

function [zobj]=DTLZ5(x)
sum3=0;
    for i=nobj-1+1:nobj-1+kk
        sum3=sum3+(x(i)-0.5)^2;
    end
    m=0;
    for n=nobj:-1:1
        m=m+1;
        zobj(m)=1;
        if n~=1
            for i=1:n-1
                zobj(m)=zobj(m)*cos((1+2*sum3*x(i))/(2*(1+sum3))*pi/2);
            end
            if n~=nobj
               zobj(m)=zobj(m)*sin((1+2*sum3*x(n))/(2*(1+sum3))*pi/2);
            else
                continue
            end
        else
            zobj(m)=zobj(m)*sin((1+2*sum3*x(n))/(2*(1+sum3))*pi/2);
        end
        zobj(m)=zobj(m)*(1+sum3);
    end
end

% DTLZ6

function [zobj]=DTLZ6(x)
sum3=0;
    for i=nobj-1+1:nobj-1+kk
        sum3=sum3+(x(i))^0.1;
    end
    m=0;
    for n=nobj:-1:1
        m=m+1;
        zobj(m)=1;
        if n~=1
            for i=1:n-1
                zobj(m)=zobj(m)*cos((1+2*sum3*x(i))/(2*(1+sum3))*pi/2);
            end
            if n~=nobj
               zobj(m)=zobj(m)*sin((1+2*sum3*x(n))/(2*(1+sum3))*pi/2);
            else
                continue
            end
        else
            zobj(m)=zobj(m)*sin((1+2*sum3*x(n))/(2*(1+sum3))*pi/2);
        end
        zobj(m)=zobj(m)*(1+sum3);
    end
end
end