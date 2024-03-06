function prob = probabilityYgiventransitions(yn,x,xprev,s,sprev,M,q01,q10)
    k = s - sprev;
    sbarprev = M-1-sprev;
    q00 = 1 - q01;
    q11 = 1 - q10;
    yplus = 0;
    yminus = 0;
    y1 = 0;
    y0 = 0;
    yI = 0;
    switch yn
        case{'I'}
            if abs(k)
                prob = 0;
            else
                if x ~= xprev
                    prob = 0;
                else
                    prob = q00^sprev*q11^(M-1-sprev)/(calculate_transition_probability(0, sprev, M, q01, q10));
                end
            end
        case{'+'}
            if x~= xprev
                prob = 0;
            else
                if (k) ~= -1
                    prob = 0;
                else
                    prob = sprev*q01*q00^(sprev-1)*q11^(sbarprev)/(calculate_transition_probability(1, sprev, M, q01, q10));
                end
            end
        case{'-'}
            if x~= xprev
                prob = 0;
            else
                if (k) ~= 1
                    prob = 0;
                else
                    prob = sbarprev*q10*q11^(sbarprev-1)*q00^sprev/(calculate_transition_probability(k, sprev, M, q01, q10));
                end
            end
        case{'1'}
            if(x~=1 || xprev ~=0 || s~=sprev)
                prob = 0;
            else
                prob = q00^sprev*q11^(sbarprev)/(calculate_transition_probability(k, sprev, M, q01, q10));
            end
        case{'0'}
            if(x~=0 || xprev ~=1 || s~=sprev)
                prob = 0;
            else
                prob = q00^sprev*q11^(sbarprev)/(calculate_transition_probability(k, sprev, M, q01, q10));
            end
        case{'C'}
            if abs(k)
                yI = 0;
            else
                if x ~= xprev
                    yI = 0;
                else
                    yI = q00^sprev*q11^(M-1-sprev)/(calculate_transition_probability(k, sprev, M, q01, q10));
                end
            end
            if x~= xprev
                yplus = 0;
            else
                if (k) ~= -1
                    yplus = 0;
                else
                    yplus = sprev*q01*q00^(sprev-1)*q11^sbarprev/(calculate_transition_probability(k, sprev, M, q01, q10));
                end
            end
            if x~= xprev
                yminus = 0;
            else
                if (k) ~= 1
                    yminus = 0;
                else
                    yminus = sbarprev*q10*q11^(sbarprev-1)*q00^sprev/(calculate_transition_probability(k, sprev, M, q01, q10));
                end
            end
            if(x~=1 || xprev ~=0 || s~=sprev)
                y1 = 0;
            else
                y1 = q00^sprev*q11^(sbarprev)/(calculate_transition_probability(0, sprev, M, q01, q10));
            end
            if(x~=0 || xprev ~=1 || s~=sprev)
                y0 = 0;
            else
                y0 = q00^sprev*q11^(sbarprev)/(calculate_transition_probability(0, sprev, M, q01, q10));
            end
            prob = 1 - (yI + yplus + yminus + y1 + y0);
    end

end