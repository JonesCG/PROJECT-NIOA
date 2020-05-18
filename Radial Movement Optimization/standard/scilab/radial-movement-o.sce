
clear
clc

exec('objective.sce');

POPSIZE = 300
DIM = 100
UPPER = 100
LOWER = -100
MAXITER = 1000
k = 1
C1 = 0.8
C2 = 0.9
WMIN = 0
WMAX = 1
VMAX = (UPPER-LOWER)/k;
disp('RUNNING...')
rand('seed',getdate('s'))

P = rand(POPSIZE,DIM).*(UPPER-LOWER)+LOWER
FIT = F1(P)
[RBESTFIT IND] = min(FIT)
RBESTPOS = P(IND,:)
GBESTPOS = P(IND,:)
GBESTFIT = RBESTFIT

CENTER = P(IND,:)
/*
xtitle("INITIALIZATION")
square(LOWER,LOWER,UPPER,UPPER)
plot(P(:,1),P(:,2),'b.')
plot(GBESTPOS(1,1),GBESTPOS(1,2),'g.')
//xs2png(gcf(),'gif/rmo/ITER0.png')
clf()
*/

disp(GBESTFIT)
for ITER=1:MAXITER
    //disp(ITER)
    VEL = grand(POPSIZE,DIM,'unf',-1,1).*VMAX
    //VEL = rand(POPSIZE,DIM)
    WEIGHT = WMAX - (WMAX-WMIN)*((ITER-1)/MAXITER)
    P = VEL.*WEIGHT + repmat(CENTER,[POPSIZE 1])
    P = ammend(P,UPPER,LOWER)
    FIT = F1(P)
    [RBESTFIT IND] = min(FIT)
    RBESTPOS = P(IND,:)
    UPDATE = (GBESTPOS-CENTER).*C1 + (RBESTPOS-CENTER).*C2
    CENTER = CENTER + UPDATE
    if RBESTFIT < GBESTFIT
        GBESTFIT = RBESTFIT
        GBESTPOS = RBESTPOS
    end
    
    FITRUN(ITER) = GBESTFIT
    /*
    xtitle("ITER "+string(ITER))
    square(LOWER,LOWER,UPPER,UPPER)
    plot(P(:,1),P(:,2),'b.')
    plot(GBESTPOS(1,1),GBESTPOS(1,2),'g.')
    //xs2png(gcf(),'gif/rmo/ITER'+string(ITER)+'.png')
    clf()
    */
end
//plot((1:MAXITER)',FITRUN,'g-')

disp(GBESTFIT)
disp('DONE')
