
clear
clc

exec('objective.sce');

rand('seed',getdate('s'))

POPSIZE = 50
DIM = 2
UPPER = 100
LOWER = -100
MAXITER = 50
KMAXITER = 50
K = 3
LAMBDA = 10

//YSGF = zeros(POPSIZE,DIM)
YSGF = rand(POPSIZE,DIM).*(UPPER-LOWER) + LOWER
FIT = F1(YSGF)

[BESTFIT IND] = min(FIT)
GBESTPOS = YSGF(IND,:)
GBESTFIT = BESTFIT

xtitle("INITIALIZATION")
square(LOWER,LOWER,UPPER,UPPER)
plot(YSGF(:,1),YSGF(:,2),'b.')
plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
xs2png(gcf(),'gif/ysga/ITER0.png')
clf()

// K-MEANS CLUSTERING
// INITIALIZE CENTROID
CENTROID = zeros(K,DIM)
CENTROID = rand(K,DIM).*(UPPER-LOWER) + LOWER

for KITER=1:KMAXITER
    POPMAT = matrix(YSGF',[1 DIM POPSIZE])
    CALC = []
    for CNUM=1:K
        CREP = repmat(CENTROID(CNUM,:),[POPSIZE 1])
        DIFF = CREP-YSGF
        DIST = sqrt(sum(DIFF.^2,'c'))
        CALC = [CALC DIST]
    end
    [VAL IND] = min(CALC,'c')
    IND_M = matrix(IND,[1 POPSIZE 1])
    VAL_M = matrix(VAL,[1 POPSIZE 1])
    for CNUM=1:K
        IND = find(IND_M == CNUM)
        if length(IND)>0
            CENTROID(CNUM,:) = mean(YSGF(IND,:),'r')
        else
            CENTROID(CNUM,:) = rand(1,DIM).*(UPPER-LOWER) + LOWER
        end
        // storing each search agent's cluster where it belongs
        ASSIGN = IND_M
    end
end
// identify chaser fish, the rest are blockers, get chaser fish index
CHASERIND = zeros(K,1)
// for each cluster
for CNUM=1:K
    // find all fish index belonging to cluster CNUM
    INDK = find(ASSIGN==CNUM)
    // get index of the fittest fish
    [VAL INDF] = min(FIT(INDK))
    // save index to chaser fish array
    CHASERIND(CNUM) = INDK(INDF)
end
//CHANGE ZONE COUNTER
Q = zeros(K,1)
// ITERATION PHASE
for ITER=1:MAXITER
    // update a
    // these parameters are for the logarithmic spiral
    a = -1-(ITER-1)*((-1-(-2))/MAXITER)
    b = 1
    // update VARIABLES for LEVY FLIGHT STEP
    ALPHA = 1
    BETA = 1.99 + (0.001)*(ITER/(MAXITER/10))
    // NUMERATOR
    NUM = gamma(1+BETA)*sin(%pi*(BETA/2))
    // DENOMINATOR
    DEN = gamma((1+BETA)/2)*BETA*2^((BETA-1)/2)
    // STANDARD DEVIATIOIN
    SIGMA_U = (NUM/DEN)^(1/BETA)
    // FOR EACH CLUSTER
    for CNUM=1:K
        // HUNTING ROUTINE FOR CHASER FISH
        // LEVY FUNCTION
        u = grand(1,DIM,"nor",0,SIGMA_U^2)
        v = grand(1,DIM,"nor",0,1)
        LEVYSTEP = u./(abs(v).^(1/BETA))
        YSGF(CHASERIND(CNUM),:) = YSGF(CHASERIND(CNUM),:) + ALPHA.*LEVYSTEP.*(YSGF(CHASERIND(CNUM),:) - GBESTPOS)
        if FIT(CHASERIND(CNUM))==GBESTPOS
            // chaser fish is the GBEST
            YSGF(CHASERIND(CNUM),:) = YSGF(CHASERIND(CNUM),:) + ALPHA.*LEVYSTEP
        end
        // BLOCKING ROUTING FOR VLOCKER FISHES
        // get indexes of blockers belonging to cluster CNUM
        INDBLOCKERS = setdiff(find(ASSIGN==CNUM),CHASERIND(CNUM))
        D = abs( rand(length(INDBLOCKERS),DIM).*repmat(YSGF(CHASERIND(CNUM),:),[length(INDBLOCKERS) 1]) -  YSGF(INDBLOCKERS,:))
        RAND = grand(length(INDBLOCKERS),DIM,"unf",a,1)
        YSGF(INDBLOCKERS,:) = D.*exp(RAND.*b).*cos(RAND.*(2*%pi)) + repmat(YSGF(CHASERIND(CNUM),:),[length(INDBLOCKERS) 1])
        YSGF(find(ASSIGN==CNUM),:) = ammend(YSGF(find(ASSIGN==CNUM),:),UPPER,LOWER)
        FIT(find(ASSIGN==CNUM)) = F1(YSGF(find(ASSIGN==CNUM),:))
        // GET FITTEST OF THE BLOCKERS
        [BESTFIT INDB] = min(FIT(INDBLOCKERS))
        if BESTFIT<FIT(CHASERIND(CNUM))
            YSGF0 = YSGF(INDBLOCKERS(INDB),:)
            YSGF(INDBLOCKERS(INDB),:) = YSGF(CHASERIND(CNUM),:)
            YSGF(CHASERIND(CNUM),:) = YSGF0
            
            FIT0 = FIT(INDBLOCKERS(INDB))
            FIT(INDBLOCKERS(INDB)) = FIT(CHASERIND(CNUM))
            FIT(CHASERIND(CNUM)) = FIT0
        end
        // IF CHASER FISH IS BETTER THAN GBEST, UPDATE
        if FIT(CHASERIND(CNUM))<GBESTFIT
            GBESTFIT = FIT(CHASERIND(CNUM))
            GBESTPOS = YSGF(CHASERIND(CNUM),:,:)
        else
            Q(CNUM) = Q(CNUM)+1
        end
        // CHANGING ZONE
        if Q(CNUM)>LAMBDA
            YSGF(find(ASSIGN==CNUM),:) = (repmat(GBESTPOS,[length(find(ASSIGN==CNUM)) 1])+YSGF(find(ASSIGN==CNUM),:))./2
            Q(CNUM) = 0
            YSGF = ammend(YSGF,UPPER,LOWER)
        end
    end
    FITRUN(ITER) = GBESTFIT
    
    //gcf().axes_size = [500 500]
    xtitle("ITER "+string(ITER))
    square(LOWER,LOWER,UPPER,UPPER)
    plot(YSGF(:,1),YSGF(:,2),'b.')
    plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
    xs2png(gcf(),'gif/ysga/ITER'+string(ITER)+'.png')
    clf()
    
end
//plot((1:MAXITER)',FITRUN,'g-')
