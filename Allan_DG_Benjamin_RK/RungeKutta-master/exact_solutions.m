function [U,flag] = exact_solutions(t,options)

flag = 0;

file_name = [options.model.function '_ex'];

if (exist(file_name,'file'))
    func = str2func(file_name);
    U    = func(t,options);
    flag = 1;
else
    U = 0;
    disp('no exact solution');
end

end
