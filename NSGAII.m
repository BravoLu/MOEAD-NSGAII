%NSGA-II
function NSGAII(Problem,M)
clc;format compact;tic;

%-----------------------------------------------------------------------------------------
%�����趨
Generations = 700;
if M == 2
    N = 100;
else M == 3
    N = 105;
end
%-----------------------------------------------------------------------------------------

%�㷨��ʼ
    %��ʼ����Ⱥ
    [Population,Boundary,Coding] = P_objective('init',Problem,M,N);
    FunctionValue = P_objective('value',Problem,M,Population);
    %This is for test -----Start
    Seq = zeros(size(FunctionValue));
    for i = 1 : M
        [~,tmp] = sort(FunctionValue(:,i));
        [~,Seq(:,i)] = sort(tmp);
    end
    %this is for test -----Start
    FrontValue = P_sort(FunctionValue);
    CrowdDistance = F_distance(FunctionValue,FrontValue);
    
    %��ʼ����
    for Gene = 1 : Generations    
        %�����Ӵ�
        MatingPool = F_mating(Population,FrontValue,CrowdDistance);
        Offspring = P_generator(MatingPool,Boundary,Coding,N);
        Population = [Population;Offspring];
        FunctionValue = P_objective('value',Problem,M,Population);
        [FrontValue,MaxFront] = P_sort(FunctionValue,'half');
        CrowdDistance = F_distance(FunctionValue,FrontValue);

        
        %ѡ����֧��ĸ���        
        Next = zeros(1,N);
        NoN = numel(FrontValue,FrontValue<MaxFront);
        Next(1:NoN) = find(FrontValue<MaxFront);
        
        %ѡ�����һ����ĸ���
        Last = find(FrontValue==MaxFront);
        [~,Rank] = sort(CrowdDistance(Last),'descend');
        Next(NoN+1:N) = Last(Rank(1:N-NoN));
        
        %��һ����Ⱥ
        Population = Population(Next,:);
        FrontValue = FrontValue(Next);
        CrowdDistance = CrowdDistance(Next);
        
		FunctionValue = P_objective('value',Problem,M,Population);
		cla;
		for i = 1 : MaxFront
			FrontCurrent = find(FrontValue==i);
			P_draw(FunctionValue(FrontCurrent,:));
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

%���ɽ��
    %P_output(Population,toc,'NSGA-II',Problem,M,Run);
end
