p01 = 0.1;
p10 = 0.01
P = [1-p01 p01;p10 1 - p10];
B = emitprobCalc(5);
lambdarange = -10:0.1:10;

lambdamat = zeros(3,length(lambdarange));

for output = 1:3
    for lambdaswiper = 1:length(lambdarange)
        lambdamat(output,lambdaswiper) = LLR_random(output,lambdarange(lambdaswiper),P,B);
    end
end

function emitprob = emitprobCalc(M)
    emitprob = [1/M*(1-1/M)^(M-1) 0 1-1/M*(1-1/M)^(M-1);0 1/M*(1-1/M)^(M-1) 1-1/M*(1-1/M)^(M-1)];
end
function lambda = LLR_random(y,prevlambda,P,B)
    lambda = log(B(1,y)) - log(B(2,y)) + log(P(1,1) + P(2,1)*exp(-prevlambda)) - log(P(2,2) + P(1,2)*exp(-prevlambda));
end
