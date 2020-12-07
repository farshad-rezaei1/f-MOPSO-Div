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
function [indexx]=Archive(indexx,nobj,znp,np_archive)
% Preserving a certain solution in the external archive once detected to be
% a non-dominated solution and removing the ones previously recorded in the
% archive once detected to be dominated by the current one

if np_archive~=1
    for y=1:np_archive-1
        if indexx(1,y)==1
            c1=0;c2=0;
            for i=1:nobj
                if (znp(i,np_archive)<=znp(i,y)) 
                    c1=c1+1;
                elseif (znp(i,np_archive)>=znp(i,y))
                    c2=c2+1;
                end
            end
            if c1==nobj 
                indexx(1,y)=0;
            elseif c2==nobj && c2~=c1
                indexx(1,np_archive)=0;
                break
            end
        end
    end
end