function [refinedProjections, thetasestimated, shiftsestimated, classestimated] = ...
    SHARPord_cluster(Pgiven, svector, sigmaNoise, shift_amplitude, initialshiftestimate,...
        initialThetaEstimate, noisyOrientations, initialClassEstimated, angle_amplitude)

    addpath(genpath('../utilities'));

    numkeep = size(Pgiven, 2);
    kmax = numkeep - 1;
    if noisyOrientations == 1
        numstarts = 2;
    else
        numstarts = 12;
    end
    Ord = 8;

    Pgiven = denoise(Pgiven, sigmaNoise, 50, 700);
    Pgiven = max(0, Pgiven);
    refinedProjections = Pgiven;

    numiter = 40;   

    thetasestimated_bystart = zeros(numkeep, numstarts);
    shiftsestimated_bystart = zeros(numkeep, numstarts);
    classestimated_bystart = zeros(numkeep, numstarts);
    Ehlccvalues_bystart = zeros(numstarts, 1);

    if noisyOrientations == 1
        min_limit = initialThetaEstimate - angle_amplitude;
        max_limit = initialThetaEstimate + angle_amplitude;
    else
        min_limit = repmat(-179, numkeep, 1);
        max_limit = repmat(180, numkeep, 1);
    end

    parfor start = 1:numstarts
        shiftestimated = initialshiftestimate;
        shiftedPgiven = correct_projection_shifts(Pgiven, shiftestimated);

        classestimated = initialClassEstimated; 
        
        if noisyOrientations == 1
            thetasestimated = initialThetaEstimate;
        else
            % Initialize angles randomly
            thetasestimated = randi([1 179], numkeep, 1);
        end

        Ehlcc = calc_error(...
            shiftedPgiven, kmax, svector, Ord, thetasestimated, classestimated);

        Ehlccvalues = zeros(numiter * (numkeep - 1), 1);
        deltas = zeros(numiter, 1);
        ctr = 1;

        Ehlccvalues(ctr) = Ehlcc;

        for iteration = 1:numiter
            disp([num2str(start), ', ', num2str(iteration)]);

            for i = 1:numkeep
                ctr = ctr + 1;

                besttheta = thetasestimated(i);
                finalbesttheta2 = thetasestimated(i);
                finalbesttheta = thetasestimated(i);
                bestshift = shiftestimated(i);
                finalbestshift = shiftestimated(i);
                bestclass = classestimated(i);

                bestshiftE = Ehlcc;
                bestE = Ehlcc;
                bestclassE = Ehlcc;
                
                for c = 1:3
                    class_iter = classestimated;
                    class_iter(i) = c;

                    for s = -shift_amplitude:shift_amplitude
                        shift_iter = shiftestimated;
                        shift_iter(i) = s;
                        
                        shiftedPgiven = ....
                            correct_projection_shifts(Pgiven, shift_iter);

                        if noisyOrientations ~= 1
                            if c == 1
                                max_limit_i = max_limit(i);
                                min_limit_i = min_limit(i);
                            else
                                max_limit_i = besttheta + 2;
                                min_limit_i = besttheta - 2;
                            end
                        else
                            max_limit_i = max_limit(i);
                            min_limit_i = min_limit(i);
                        end
                         
                        err_t = err_for_all_angles(...
                            shiftedPgiven, kmax, svector, Ord, thetasestimated, class_iter,...
                            max_limit_i, min_limit_i, i);
                        
                        [E_t, idx_err_t] = min(err_t);
                        if E_t < bestE
                            bestE = E_t;
                            besttheta = idx_err_t + min_limit(i) - 1;
                        end
                        
                        if bestE < bestshiftE
                            bestshift = -s;
                            bestshiftE = bestE;
                            finalbesttheta = besttheta;
                        end
                    end

                    if bestshiftE < bestclassE
                        bestclassE = bestshiftE;
                        bestclass = c;
                        finalbestshift = bestshift;
                        finalbesttheta2 = finalbesttheta;
                    end
                end

                thetasestimated(i) = finalbesttheta2;
                shiftestimated(i) = finalbestshift;
                classestimated(i) = bestclass;
                
                shiftedPgiven = correct_projection_shifts(Pgiven, shiftestimated);
                Ehlcc = calc_error(shiftedPgiven, kmax, svector, Ord, thetasestimated, classestimated);
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
        classestimated_bystart(:, start) = classestimated;
    end % multistart loop

    [~, optstart] = min(Ehlccvalues_bystart);
    thetasestimated = thetasestimated_bystart(:, optstart);
    shiftsestimated = shiftsestimated_bystart(:, optstart);
    classestimated = classestimated_bystart(:, optstart);
end
