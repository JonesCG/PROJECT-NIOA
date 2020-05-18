
clear
clc

exec('objective.sce');

disp("RUNNING...")
rand('seed',getdate('s'))

POPSIZE = 300
DIM = 100
UPPER = 100
LOWER = -100
MAXITER = 1000
WEIGHT = 1
C1 = 0.8
C2 = 0.9

P = rand(POPSIZE,DIM).*(UPPER-LOWER) + LOWER
VEL = zeros(POPSIZE,DIM)
FIT = F1(P)

[BESTFIT IND] = min(FIT)
PBESTPOS = P
PBESTFIT = FIT
GBESTPOS = P(IND,:)
GBESTFIT = FIT(IND)
INITBEST = GBESTFIT

/*
xtitle("INITIALIZATION")
square(LOWER,LOWER,UPPER,UPPER)
plot(P(:,1),P(:,2),'b.')
plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
xs2png(gcf(),'gif/pso/ITER0.png')
clf()
*/

for ITER=1:MAXITER
    //disp(ITER)
    WEIGHT = 1 - (1-0)*((ITER-1)/MAXITER)
    PB = rand(POPSIZE,DIM).*(PBESTPOS-P)
    GB = rand(POPSIZE,DIM).*(repmat(GBESTPOS,[POPSIZE 1])-P)
    VEL = VEL.*WEIGHT + PB.*C1 + GB.*C2
    P = P + VEL
    P = ammend(P,UPPER,LOWER)
    FIT = F1(P)
    I = find(FIT<PBESTFIT)
    PBESTPOS(I,:) = P(I,:)
    PBESTFIT(I) = FIT(I)
    [BESTFIT IND] = min(PBESTFIT)
    if BESTFIT < GBESTFIT
        GBESTFIT = BESTFIT
        GBESTPOS = PBESTPOS(IND,:)
    end
    FITRUN(ITER) = GBESTFIT
    
    /*
    xtitle("ITER "+string(ITER))
    square(LOWER,LOWER,UPPER,UPPER)
    plot(P(:,1),P(:,2),'b.')
    plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
    xs2png(gcf(),'gif/pso/ITER'+string(ITER)+'.png')
    clf()
    */
    
end
//plot((1:MAXITER)',FITRUN,'g-')
disp(INITBEST)
disp(GBESTFIT)
disp('DONE')
