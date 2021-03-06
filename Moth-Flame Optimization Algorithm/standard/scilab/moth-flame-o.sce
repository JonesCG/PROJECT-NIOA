
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

b=1

MOTH = rand(POPSIZE,DIM).*(UPPER-LOWER)+LOWER
FIT = F1(MOTH)
[BESTFIT IND] = min(FIT)
GBESTPOS = MOTH(IND,:)
GBESTFIT = BESTFIT
// SORT MOTH BASED ON FITNESS AND STORE SORTED MOTH TO FLAME
[FLAMEFIT SORTIND] = gsort(FIT,'lr','i')
FLAME = MOTH(SORTIND,:,:)

PREVPOP = MOTH
PREVFIT = FIT

xtitle("INITIALIZATION")
square(LOWER,LOWER,UPPER,UPPER)
plot(MOTH(:,1),MOTH(:,2),'b.')
plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
xs2png(gcf(),'gif/mfo/ITER0.png')
clf()
for ITER=1:MAXITER
    // VALUE OF a, decreases linearly from -1 to -2
    a = -1 + (ITER-1)*((-2-(-1))/MAXITER)
    // UPDATE FLAMENUMBER
    FLAMENUM = round(POPSIZE - ITER*((POPSIZE-1)/MAXITER))
    t = grand(POPSIZE,DIM,"unf",a,1)
    // CALCULATE DISTANCE OF THE ith MOTH TO jth FLAME
    DTF = abs(FLAME - MOTH)
    E = exp(b.*t)
    COSF = cos(t.*2*%pi)
    // UPDATE MOTHS TOWARDS TO THEIR CORRESPONDING FLAMES
    MOTH = DTF.*E.*COSF + FLAME
    // PROGRAM ENTERS HERE IF THERE ARE MOTHS THAT DOES NOT BELONG TO THE TOP 1 TO FLAMENUM
    if FLAMENUM < POPSIZE 
        // CALCULATE DISTANCE of the moth to the WORST FLAME
        DTF = abs(repmat(FLAME(FLAMENUM,:),[(POPSIZE-FLAMENUM) 1]) - MOTH((FLAMENUM+1:POPSIZE),:))
        E = exp(b.*t((FLAMENUM+1:POPSIZE),:))
        COSF = cos(t((FLAMENUM+1:POPSIZE),:).*2.*%pi)
        // UPDATE MOTHS TOWARDS THE WORST FLAME
        MOTH((FLAMENUM+1:POPSIZE),:) = DTF.*E.*COSF + FLAME((FLAMENUM+1:POPSIZE),:)
    end
    
    MOTH = ammend(MOTH,UPPER,LOWER)
    FIT = F1(MOTH)
    // COMBINE PREVPOP WITH FLAME AND PREVFIT WITH FLAMEFIT
    DOUBLEPOP = [PREVPOP;FLAME]
    DOUBLEFIT = [PREVFIT;FLAMEFIT]
    // SORT DOUBLEPOP ACCORDING TO FITNESS
    [DFSORTED SORTIND] = gsort(DOUBLEFIT,'lr','i')
    DPSORTED = DOUBLEPOP(SORTIND,:)
    // UPDATE FLAMES
    FLAMEFIT = DFSORTED(1:POPSIZE)   // GET TOP 1 to POPSIZE
    FLAME = DPSORTED(1:POPSIZE,:)   // GET TOP 1 to POPSIZE
    // UPDATE PREVIOUS POP
    PREVPOP = MOTH
    PREVFIT = FIT
    if FLAMEFIT(1)<GBESTFIT
        GBESTFIT = FLAMEFIT(1)
        GBESTPOS = FLAME(1,:)
    end
    FITRUN(ITER) = GBESTFIT
    
    //gcf().axes_size = [500 500]
    xtitle("ITER "+string(ITER))
    square(LOWER,LOWER,UPPER,UPPER)
    plot(MOTH(:,1),MOTH(:,2),'b.')
    plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
    xs2png(gcf(),'gif/mfo/ITER'+string(ITER)+'.png')
    clf()
    
end
