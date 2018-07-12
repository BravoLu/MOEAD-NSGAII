function MOEAD(Problem,M)
clc;format compact;tic;


%�����趨
Generations = 700;
delta = 0.9;
nr = 2;
if M == 2
    N = 100;
    H = 99;
else M == 3
    N = 105;
    H = 13;
end

    %��ʼ������
    Evaluations = Generations*N;
    [N,W] = EqualWeight(H,M);
    W(W==0) = 0.000001;
    T = floor(N/10);
    Generations = floor(Evaluations/N);

    %�ھ��ж�
    B = zeros(N);
    for i = 1 : N-1
        for j = i+1 : N
            B(i,j) = norm(W(i,:)-W(j,:));
            B(j,i) = B(i,j);
        end
    end
    [~,B] = sort(B,2);
    B = B(:,1:T);
    
    %��ʼ����Ⱥ
    [Population,Boundary] = Objective(0,Problem,M,N);
    FunctionValue = Objective(1,Problem,M,Population);
    Z = min(FunctionValue);

    %��ʼ����
    for Gene = 1 : Generations
        %��ÿ������ִ�в���
        for i = 1 : N
            %ѡ����ĸ
            if rand < delta
                P = B(i,:);
            else
                P = 1:N;
            end
            k = randperm(length(P));
            
            %�����Ӵ�
            Offspring = Gen(Population(i,:),Population(P(k(1)),:),Population(P(k(2)),:),Boundary);
            OffFunValue = Objective(1,Problem,M,Offspring);

            %�������������
            Z = min(Z,OffFunValue);
            
            %����P�еĸ���
            c = 0;
            for j = randperm(length(P))
                if c >= nr
                    break;
                end
                g_old = max(abs(FunctionValue(P(j),:)-Z).*W(P(j),:));
                g_new = max(abs(OffFunValue-Z).*W(P(j),:));              
                if g_new < g_old
                    %���µ�ǰ�����ĸ���
                    Population(P(j),:) = Offspring;
                    FunctionValue(P(j),:) = OffFunValue;
                    c = c+1;
                end
            end

        end
        cla;
        DrawGraph(FunctionValue);
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
        %clc;
    end
end