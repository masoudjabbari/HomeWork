function [fval] = bertrand(p)
    % Set parameters:
    v = [2; 2];
%     these parameters should be incoming parameters in a function
    
%     end every line with semicolumn, otherwise it prints everything in
%     command line. only the answer to a question should be printed. 
   exp_excess=exp(v-p)
   D=exp_excess./(1+sum(exp_excess))

    %Evaluate function: 
    fval = ones(size(p))-p+p.*D
end