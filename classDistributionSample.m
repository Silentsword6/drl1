% Date Created: 8/25/10
% Description:
% This is a distribution sampler designed to pick random variable values
% from the domain of a probabilty mass function (pmf) by comparing to the
% cummulative distribution function for a continuous uniform distribution 
% from 0 to 1.
%这是一种分布采样器，通过与从0到1的连续均匀分布的累积分布函数进行比较，从概率质量函数(pmf)的域中挑选随机变量值。

function [measure] = classDistributionSample(cm, tgt)

    % importance sample to generate a measurement from the given track
    r = rand(1);
    
    total = 0;
    for m = 1:size(cm,1)
        
        total = total + cm(m,tgt);
        
        if r <= total 
            
            measure = m;
            break;
            
        end
        
    end

end
