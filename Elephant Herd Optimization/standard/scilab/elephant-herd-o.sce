
clear
clc

exec('objective.sce');

disp("RUNNING...")
rand('seed',getdate('s'))

POPSIZE = 50
DIM = 2
UPPER = 100
LOWER = -100
MAXITER = 50
// SCALE FACTOR ALPHA a and BETA b
a = 0.5
b = 0.1
// NUMBER OF CLAN
NCLAN = 5

ELEPHANT = rand(POPSIZE,DIM).*(UPPER-LOWER) + LOWER
FIT = F1(ELEPHANT)
// GET GBEST
[BESTFIT IND] = min(FIT)
GBESTFIT = BESTFIT
GBESTPOS = ELEPHANT(IND,:)

xtitle("INITIALIZATION")
square(LOWER,LOWER,UPPER,UPPER)
plot(ELEPHANT(:,1),ELEPHANT(:,2),'b.')
plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
xs2png(gcf(),'gif/eho/ITER0.png')
clf()

for ITER=1:MAXITER
    [FIT SORTIND] = gsort(FIT,'lr','i')
    ELEPHANT = ELEPHANT(SORTIND,:)
    // CLAN UPDATING OPERATOR
    estart = 1
    member = POPSIZE/NCLAN
    for c=1:NCLAN
        IN = estart:member*c
        ELEPHANT(IN,:) = ELEPHANT(IN,:) + a.*(repmat(ELEPHANT(estart,:),[member 1]) - ELEPHANT(IN,:)).*rand(member,DIM)
        CENTER = mean(ELEPHANT(estart:member*c,:),'r')
        ELEPHANT(estart,:) = CENTER.*b
        estart = estart+member
    end
    
    // SEPARATING OPERATOR
    WN = matrix(1:POPSIZE,[member NCLAN])(member,:)
    ELEPHANT(WN,:) = LOWER + (UPPER - LOWER +1.0).*rand(NCLAN,DIM)
    
    ELEPHANT = ammend(ELEPHANT,UPPER,LOWER)
    // EVALUATE
    FIT = F1(ELEPHANT)
    [BESTFIT IND] = min(FIT)
    if BESTFIT < GBESTFIT
        GBESTFIT = BESTFIT
        GBESTPOS = ELEPHANT(IND,:)
    end
    FITRUN(ITER) = GBESTFIT
    
    xtitle("ITER "+string(ITER))
    square(LOWER,LOWER,UPPER,UPPER)
    plot(ELEPHANT(:,1),ELEPHANT(:,2),'b.')
    plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
    xs2png(gcf(),'gif/eho/ITER'+string(ITER)+'.png')
    clf()
    
end
//plot((1:MAXITER)',FITRUN,'g-')
