function [refinedProjections, thetasestimated, shiftsestimated] = ...
    SHARPord(Pgiven, svector, sigmaNoise, shift_amplitude, initialshiftestimate,...
        initialThetaEstimate, noisyOrientations)

    addpath(genpath('../utilities'));

    numkeep = size(Pgiven, 2);
    kmax = numkeep - 1;
    if noisyOrientations == 1
        numstarts = 2;
    else
        numstarts = 20;
    end
    Ord = 8;

    Pgiven = denoise(Pgiven, sigmaNoise, 50, 700);
    Pgiven = max(0, Pgiven);
    refinedProjections = Pgiven;

    numiter = 40;   

    thetasestimated_bystart = zeros(numkeep, numstarts);
    shiftsestimated_bystart = zeros(numkeep, numstarts);
    Ehlccvalues_bystart = zeros(numstarts, 1);

    if noisyOrientations == 1
        min_limit = initialThetaEstimate - 1;
        max_limit = initialThetaEstimate + 1;
    else
        min_limit = repmat(-179, numkeep, 1);
        max_limit = repmat(180, numkeep, 1);
    end

    parfor start = 1:numstarts
        shiftestimated = initialshiftestimate;
        shiftedPgiven = correct_projection_shifts(Pgiven, shiftestimated);

        PMord = ...
            assemblePMord(shiftedPgiven, kmax, svector, Ord, numkeep);
        
        if noisyOrientations == 1
            thetasestimated = initialThetaEstimate;
        else
            % Initialize angles randomly
            thetasestimated = randi([1 179], numkeep, 1);
        end

        A = assembleA(thetasestimated, Ord); 
        IMestimated = A \ PMord;

        Ehlccvec = A * IMestimated - PMord; 
        Ehlcc = norm(Ehlccvec, 2);
        Ehlccvalues = zeros(numiter * (numkeep - 1), 1);
        deltas = zeros(numiter, 1);
        ctr = 1;

        Ehlccvalues(ctr) = Ehlcc;

        for iteration = 1:numiter
            disp([num2str(start), ', ', num2str(iteration)]);

            for i = 1:numkeep
                ctr = ctr + 1;

                besttheta = thetasestimated(i);
                finalbesttheta = thetasestimated(i);
                bestshift = shiftestimated(i);
                bestshiftE = Ehlcc;
                bestE = Ehlcc;
                
                for s = -shift_amplitude:shift_amplitude
                    shift_iter = shiftestimated;
                    shift_iter(i) = s;
                    
                    shiftedPgiven = ....
                        correct_projection_shifts(Pgiven, shift_iter);
                    PMord = ...
                        assemblePMord(shiftedPgiven, kmax, svector, Ord, numkeep);
                    for t = min_limit(i):max_limit(i)

                        thetas_iter = thetasestimated;
                        thetas_iter(i) = t;
                        A = assembleA(thetas_iter, Ord);
                        IMestimated = A \ PMord;
                        E_tvec = A * IMestimated - PMord;
                        E_t = norm(E_tvec, 2);

                        if E_t < bestE
                            besttheta = t;
                            bestE = E_t;
                        end
                    end
                    
                    if bestE < bestshiftE
                        bestshift = -s;
                        bestshiftE = bestE;
                        finalbesttheta = besttheta;
                    end
                end

                thetasestimated(i) = finalbesttheta;
                shiftestimated(i) = bestshift;
                
                shiftedPgiven = correct_projection_shifts(Pgiven, shiftestimated);
                PMord = ...
                    assemblePMord(shiftedPgiven, kmax, svector, Ord, numkeep);

                A = assembleA(thetasestimated, Ord);
                IMestimated = A \ PMord;

                Ehlccvec = A * IMestimated - PMord;
                Ehlcc = norm(Ehlccvec, 2);
                Ehlccvalues(ctr) = Ehlcc;

            end
            Ehlcciter = Ehlccvalues(ctr);
            Ehlccpreviter = Ehlccvalues(ctr - (numkeep)); % -1 iff theta0 is not updated

            delta = Ehlccpreviter - Ehlcciter;
            deltas(iteration) = delta;
            if abs(delta) < 0.1
                break
            end

        end
        Ehlccvalues_bystart(start) = Ehlccvalues(ctr);
        thetasestimated_bystart(:, start) = thetasestimated;
        shiftsestimated_bystart(:, start) = shiftestimated;
    end % multistart loop

    [~, optstart] = min(Ehlccvalues_bystart);
    thetasestimated = thetasestimated_bystart(:, optstart);
    shiftsestimated = shiftsestimated_bystart(:, optstart);
end
