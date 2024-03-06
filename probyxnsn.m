function logprob = probyxnsn(y,n,xn,sn,M,q01,q10)
    global Y_history
    global alphaTable
    yprevIndex = mapYtoIndex(Y_history(n-1));
    q00 = 1-q01;
    q11 = 1 - q10;
    prob = 0;
    logprob = -Inf;
    for xprev = 0:1
        for sprev = 0:M-1
            %compute P(Yn|\simga_{n} and \sigma_{n-1})
            P1 = probabilityYgiventransitions(y,xn,xprev,sn,sprev,M,q01,q10);
            if P1 == 0
                continue
            end
            %compute P(\sigma_{n} | \sigma_{n-1}) = A1 * A2, where A1 = P(x_n|x_{n-1}) and A2 =  P(\sigma_n|\sigma_{n-1})
            k = sn - sprev;
            A2 = calculate_transition_probability(k, sprev, M, q01, q10);
            if xn == 0
                xstring = '0';
            else
                xstring = '1';
            end            
            if xprev == 0
                xstring = strcat(xstring, '0');
            else
                xstring = strcat(xstring, '1');
            end
            switch xstring
                case{'00'}
                    A1 = q00;
                case{'01'}
                    A1 = q10;
                case{'10'}
                    A1 = q01;
                case{'11'}
                    A1 = q11;
            end
            P2 = A1*A2;
            if P2 == 0
                continue
            end
            %compute the multiplication, given the previous term            
            %Pt = P1*P2*alphaTable(xprev+1,sprev+1,n-1);
            log_Pt = log(P1) + log(P2) + alphaTable(xprev+1,sprev+1,n-1);
            %prob = prob + Pt;
            logprob = max(logprob, log_Pt) + log(1 + exp(-abs(logprob - log_Pt)));
        end
    end

end