
clear
clc

exec('objective.sce');

POPSIZE = 300
DIM = 100
UPPER = 100
LOWER = -100
MAXITER = 1000

disp("RUNNING...")

rand('seed',getdate('s'))

STAR = rand(POPSIZE,DIM).*(UPPER-LOWER) + LOWER
FIT = F1(STAR)

[BHFIT IND] = min(FIT)
BH = STAR(IND,:)
GBESTFIT = BHFIT
GBESTPOS = BH
/*
xtitle("INITIALIZATION")
square(LOWER,LOWER,UPPER,UPPER)
plot(STAR(:,1),STAR(:,2),'b.')
//plot(BH(:,1),BH(:,2),'kx')
plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
//xs2png(gcf(),'gif/bhoa/ITER0.png')
clf()
*/

disp(GBESTFIT)

for ITER=1:MAXITER
    //disp(ITER)
    STAR = STAR + rand(POPSIZE,DIM).*(repmat(BH,[POPSIZE 1]) - STAR)
    STAR = ammend(STAR,UPPER,LOWER)
    FIT = F1(STAR)
    [BESTFIT IND] = min(FIT)
    if BESTFIT < BHFIT
        
        TMPPOS = STAR(IND,:)
        TMPFIT = FIT(IND)
        STAR(IND,:) = BH
        FIT(IND) = BHFIT
        BH = TMPPOS
        BHFIT = TMPFIT
        
        /*
        STAR(IND,:) = STAR(IND,:)+BH
        FIT(IND) = FIT(IND)+BHFIT
        BH = STAR(IND,:)-BH
        BHFIT= FIT(IND)-BHFIT
        STAR(IND,:) = STAR(IND,:)-BH
        FIT(IND) = FIT(IND)-BHFIT
        */
    end
    
    R = BHFIT/sum(FIT)
    D = sqrt( sum( (repmat(BH,[POPSIZE 1])-STAR).^2,'c') )
    INDS = find(D<R)
    if length(INDS)>0
        STAR(INDS,:,1) = rand(length(INDS),DIM).*(UPPER-LOWER)+LOWER
    end
    
    if BHFIT<GBESTFIT
        GBESTFIT = BHFIT
        GBESTPOS = BH
    end
    
    /*
    xtitle("ITER "+string(ITER))
    square(LOWER,LOWER,UPPER,UPPER)
    plot(STAR(:,1),STAR(:,2),'b.')
    //plot(BH(:,1),BH(:,2),'kx')
    plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
    //xs2png(gcf(),'gif/bhoa/ITER'+string(ITER)+'.png')
    clf()
    */
    
    FITRUN(ITER)=GBESTFIT
end
//plot((1:MAXITER)',FITRUN,'k-')
disp(GBESTFIT)
disp("DONE")
