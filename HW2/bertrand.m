function [fval] = bertrand(p)
    % Set parameters:
    v = [2; 2];
   exp_excess=exp(v-p)
   D=exp_excess./(1+sum(exp_excess))

    %Evaluate function: 
    fval = ones(size(p))-p+p.*D
end