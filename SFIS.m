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
function [dom,dom_ppbest]=SFIS(np_archive,nobj,zobj,z_ppbest,znp,z_min,z_max,st,it)
np=size(z_ppbest,2);
mf=zeros(1,nobj);
nweight=factorial(1/st+nobj-1)/factorial(nobj-1)/factorial(1/st+nobj-1-nobj+1);
sigma_z=zeros(2,nobj);
mu_z=zeros(2,nobj);
z=zeros(nobj,np);
w=zeros(nobj);
mu=zeros(1,nobj*nweight);
zz=zeros(1,nobj*nweight);
dom=zeros(np);
dom_ppbest=zeros(np);

% Deriving the Statistical Parameters (Std and Mean values of the High and Low classes of the objectives) 
for i = 1:nobj
    sigma_z(1,i)=std(znp(i,1:floor(np_archive/2)))*(1/(z_max(i)-z_min(i)));
    sigma_z(2,i)=std(znp(i,floor(np_archive/2)+1:np_archive))*(1/(z_max(i)-z_min(i)));
    mu_z(1,i)=(mean(znp(i,1:floor(np_archive/2)))-z_max(i))/(z_max(i)-z_min(i));
    mu_z(2,i)=(mean(znp(i,floor(np_archive/2)+1:np_archive))-z_max(i))/(z_max(i)-z_min(i));
end
if it==1
    i3=1;
else
    i3=2;
end
j=1;

% Calculating the Comprehensive Dominance Index for the Current Particles and their Pbests
while j<=np
    for i2=1:i3
        if i2==1
            for m=1:nobj
                z(m,j)=(zobj(m)-z_min(m))/(z_max(m)-z_min(m));
            end
        elseif i2==2
            for m=1:nobj
                z(m,j)=(z_ppbest(m,j)-z_min(m))/(z_max(m)-z_min(m));
            end
        end
        sum3=0;sum4=0;
        for jj=1:nobj
            w(jj)=0;
        end
        i=0;
        for m=1:nobj
            w(m)=0.5;
            for m1=1:m-1
                w(m1)=0;
            end
            for j1=m+1:nobj
                w(j1)=0.5;
                for j2=m+1:nobj
                    if j2~=j1
                        w(j2)=0;
                    end
                end
            
                % Calculating mfs
                
                i=i+1;
                mu(1,i) = 1;
                for ii=1:nobj
                    if w(ii)==0.5
                        mf(ii)=1/(1+exp(4/(sigma_z(2,ii)*sqrt(exp(1)))*(z(ii,j)-mu_z(2,ii))));
                        mu(1,i)=mu(1,i)*mf(ii);
                    else
                        mf(ii)=1/(1+exp(4/(sigma_z(1,ii)*sqrt(exp(1)))*(z(ii,j)-mu_z(1,ii))));
                        mu(1,i)=mu(1,i)*mf(ii);
                    end
                end
                zz(1,i)=0;
                for jjj=1:nobj
                    zz(1,i)=zz(1,i)+w(jjj)*z(jjj,j);
                end
                sum3=sum3+mu(1,i)*zz(1,i);
                sum4=sum4+mu(1,i);
                i=i+1;
                mu(1,i)=1;
                for ii=1:nobj
                    mf(ii)=1/(1+exp(4/(sigma_z(2,ii)*sqrt(exp(1)))*(z(ii,j)-mu_z(2,ii))));
                    mu(1,i)=mu(1,i)*mf(ii);
                end
                zz(1,i)=0;
                for jjj=1:nobj
                    zz(1,i)=zz(1,i)+w(jjj)*z(jjj,j);
                end
                sum3=sum3+mu(1,i)*zz(1,i);
                sum4=sum4+mu(1,i);
            end
        end
        for m=1:nobj
            w(m)=1;
            for m1=1:nobj
                if m1~=m
                    w(m1)=0;
                end
            end
            
            % Calculating mfs
        
            i=i+1;
            mu(1,i)=1;
            for ii=1:nobj
                if w(ii)==1
                    mf(ii)=1/(1+exp(4/(sigma_z(2,ii)*sqrt(exp(1)))*(z(ii,j)-mu_z(2,ii))));
                    mu(1,i)=mu(1,i)*mf(ii);
                else
                    mf(ii)=1/(1+exp(4/(sigma_z(1,ii)*sqrt(exp(1)))*(z(ii,j)-mu_z(1,ii))));
                    mu(1,i)=mu(1,i)*mf(ii);
                end
            end
            zz(1,i)=0;
            for jjj=1:nobj
                zz(1,i)=zz(1,i)+w(jjj)*z(jjj,j);
            end
            sum3=sum3+mu(1,i)*zz(1,i);
            sum4=sum4+mu(1,i);
            i=i+1;
            mu(1,i)=1;
            for ii=1:nobj
                mf(ii)=1/(1+exp(4/(sigma_z(2,ii)*sqrt(exp(1)))*(z(ii,j)-mu_z(2,ii))));
                mu(1,i)=mu(1,i)*mf(ii);
            end
            zz(1,i)=0;
            for jjj=1:nobj
                zz(1,i)=zz(1,i)+w(jjj)*z(jjj,j);
            end
            sum3=sum3+mu(1,i)*zz(1,i);
            sum4=sum4+mu(1,i);
        end
        if i2==1
            dom(j)=sum3/(sum4+eps);
        elseif i2==2
            dom_ppbest(j)=sum3/(sum4+eps);
        end
        j=j+1;
    end
end