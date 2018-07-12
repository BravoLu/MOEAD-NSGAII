%NSGA-II
function NSGAII(Problem,M)
clc;format compact;tic;

%-----------------------------------------------------------------------------------------
%参数设定
Generations = 700;
if M == 2
    N = 100;
else M == 3
    N = 105;
end
%-----------------------------------------------------------------------------------------

%算法开始
    %初始化种群
    [Population,Boundary,Coding] = Objective('init',Problem,M,N);
    FunctionValue = Objective('value',Problem,M,Population);

    FrontValue = P_sort(FunctionValue);
    CrowdDistance = CrowdDistances(FunctionValue,FrontValue);
    
    %开始迭代
    for Gene = 1 : Generations    
        %产生子代
        MatingPool = F_mating(Population,FrontValue,CrowdDistance);
        Offspring = NSGA_Gen(MatingPool,Boundary,Coding,N);
        Population = [Population;Offspring];
        FunctionValue = Objective('value',Problem,M,Population);
        [FrontValue,MaxFront] = P_sort(FunctionValue,'half');
        CrowdDistance = CrowdDistances(FunctionValue,FrontValue);

        
        %选出非支配的个体        
        Next = zeros(1,N);
        NoN = numel(FrontValue,FrontValue<MaxFront);
        Next(1:NoN) = find(FrontValue<MaxFront);
        
        %选出最后一个面的个体
        Last = find(FrontValue==MaxFront);
        [~,Rank] = sort(CrowdDistance(Last),'descend');
        Next(NoN+1:N) = Last(Rank(1:N-NoN));
        
        %下一代种群
        Population = Population(Next,:);
        FrontValue = FrontValue(Next);
        CrowdDistance = CrowdDistance(Next);
        
		FunctionValue = Objective('value',Problem,M,Population);
		cla;
		for i = 1 : MaxFront
			FrontCurrent = find(FrontValue==i);
			DrawGraph(FunctionValue(FrontCurrent,:));
			hold on;
		switch Problem
			case 'DTLZ1'
				if M == 2
            		pareto_x = linspace(0,0.5);
            		pareto_y = 0.5 - pareto_x;
            		plot(pareto_x, pareto_y, 'r');
				elseif M == 3
					[pareto_x,pareto_y]  = meshgrid(linspace(0,0.5));
					pareto_z = 0.5 - pareto_x - pareto_y;
					axis([0,1,0,1,0,1]);
					mesh(pareto_x, pareto_y, pareto_z);
				end
			otherwise
				if M == 2
					pareto_x = linspace(0,1);
					pareto_y = sqrt(1-pareto_x.^2);
					plot(pareto_x, pareto_y, 'r');
				elseif M == 3
					[pareto_x,pareto_y,pareto_z] =sphere(50);
					axis([0,1,0,1,0,1]);
					mesh(1*pareto_x,1*pareto_y,1*pareto_z);
				end
		end
		pause(0.01);
		end
        clc;
        
    end
end
