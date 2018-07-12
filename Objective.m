function [Output,Boundary,Coding] = Objective(Operation,Problem,M,Input)
    persistent K;
    Boundary = NaN; Coding = NaN;
    switch Operation
        case 'init'
            k = find(~isstrprop(Problem,'digit'),1,'last');
            K = [5 10 10 10];
            K = K(str2double(Problem(k+1:end)));
            
            D = M+K-1;
            MaxValue   = ones(1,D);
            MinValue   = zeros(1,D);
            Population = rand(Input,D);
            Population = Population.*repmat(MaxValue,Input,1)+(1-Population).*repmat(MinValue,Input,1);
            
            Output   = Population;
            Boundary = [MaxValue;MinValue];
            Coding   = 'Real';
        case 'value'
            Population    = Input;
            FunctionValue = zeros(size(Population,1),M);
            switch Problem
                case 'DTLZ1'
                    g = 100*(K+sum((Population(:,M:end)-0.5).^2-cos(20.*pi.*(Population(:,M:end)-0.5)),2));
                    for i = 1 : M
                        FunctionValue(:,i) = 0.5.*prod(Population(:,1:M-i),2).*(1+g);
                        if i > 1
                            FunctionValue(:,i) = FunctionValue(:,i).*(1-Population(:,M-i+1));
                        end
                    end
                case 'DTLZ2'
                    g = sum((Population(:,M:end)-0.5).^2,2);
                    for i = 1 : M
                        FunctionValue(:,i) = (1+g).*prod(cos(0.5.*pi.*Population(:,1:M-i)),2);
                        if i > 1
                            FunctionValue(:,i) = FunctionValue(:,i).*sin(0.5.*pi.*Population(:,M-i+1));
                        end
                    end
                case 'DTLZ3'
                    g = 100*(K+sum((Population(:,M:end)-0.5).^2-cos(20.*pi.*(Population(:,M:end)-0.5)),2));
                    for i = 1 : M
                        FunctionValue(:,i) = (1+g).*prod(cos(0.5.*pi.*Population(:,1:M-i)),2);
                        if i > 1
                            FunctionValue(:,i) = FunctionValue(:,i).*sin(0.5.*pi.*Population(:,M-i+1));
                        end
                    end
                case 'DTLZ4'
                    Population(:,1:M-1) = Population(:,1:M-1).^100;
                    g = sum((Population(:,M:end)-0.5).^2,2);
                    for i = 1 : M
                        FunctionValue(:,i) = (1+g).*prod(cos(0.5.*pi.*Population(:,1:M-i)),2);
                        if i > 1
                            FunctionValue(:,i) = FunctionValue(:,i).*sin(0.5.*pi.*Population(:,M-i+1));
                        end
                    end
            end
            Output = FunctionValue;
        case 'true'
            switch Problem
                case 'DTLZ1'
                    Population = T_uniform(Input,M);
                    Population = Population/2;
                case {'DTLZ2','DTLZ3','DTLZ4'}
                Population = T_uniform(Input,M);
                    for i = 1 : size(Population,1)
                        k    = find(Population(i,:)~=0,1);
                        Temp = Population(i,[1:k-1,k+1:end])./Population(i,k);
                        Population(i,k) = sqrt(1/(sum(Temp.^2)+1));
                        Population(i,[1:k-1,k+1:end]) = Temp*Population(i,k);
                    end
            end
            Output = Population;
    end
end