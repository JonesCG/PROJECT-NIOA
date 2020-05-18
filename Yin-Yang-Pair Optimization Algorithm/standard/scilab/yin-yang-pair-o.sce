
clear
clc

exec('objective.sce');

disp("RUNNING...")
rand('seed',getdate('s'))

//POPSIZE = 2
DIM = 2
UPPER = 100
LOWER = -100
MAXITER = 50

ARCHMIN = 2
ARCHMAX = 4
//ARCHLIMIT = grand(1,1,"uin",ARCHMIN,ARCHMAX)
ARCHLIMIT = floor(rand()*(ARCHMAX+1-ARCHMIN))+ARCHMIN
ARCHLIMIT = 5
ALPHA = 10
DELTA(1:2) = 0.5

P = rand(2,DIM).*(UPPER-LOWER) + LOWER
FIT = F1(P)
[BESTFIT IND] = min(FIT)
GBESTFIT = BESTFIT
GBESTPOS = P(IND,:)
//P = normalize(P,UPPER,LOWER)
P = (P-LOWER)./(UPPER-LOWER)
ARCHBEST = []
ARCHFIT = []
ARCHCOUNT = 0

xtitle("INITIALIZATION")
square(LOWER,LOWER,UPPER,UPPER)
plot(P(:,1),P(:,2),'b.')
plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
xs2png(gcf(),'gif/yypo/ITER0.png')
clf()

for ITER=1:MAXITER
    //disp(ITER)
    // if P2 is better than P1
    if FIT(2)<FIT(1)
        // swap P1 and P2
        P0 = P(1,:)
        P(1,:) = P(2,:)
        P(2,:) = P0
        // swap their delta also
        D0 = DELTA(1)
        DELTA(1) = DELTA(2)
        DELTA(2) = D0
        // swap their fitness
        FIT0 = FIT(1)
        FIT(1) = FIT(2)
        FIT(2) = FIT0
    end
    
    // add to archive P1 and P2, their fitness, update archive count
    ARCHBEST = [ARCHBEST; P]
    ARCHFIT = [ARCHFIT; FIT]
    ARCHCOUNT = ARCHCOUNT+1
    col = ['c.' 'm.']
    for point=1:2
        //PP = actualize(P(point,:,:),ubx,uby,lbx,lby)
        //plot(PP(1,:,1),PP(1,:,2),COL(point))
        if rand()<0.5
            // 1-WAY split
            // eq. 1
            // generate 2D or 2p copies of Ppoint
            // C1 is for eq.1.1, p copies
            C1 = repmat(P(point,:),[DIM 1]) + rand(DIM,DIM).*DELTA(point)
            // C1 is for eq.1.2, p copies
            C2 = repmat(P(point,:),[DIM 1]) - rand(DIM,DIM).*DELTA(point)
            // combine generated copies, 2p copies all
            S = [C1; C2]
        else
            // create the random binary matrix - no bitstrings are equal
            VALPRM = grand(1,"prm",(1:2^DIM))
            // Get the first 2D value
            DVAL = VALPRM(1:2*DIM)
            // convert to binary (string)
            BINSTR = dec2bin(DVAL-1,DIM)
            // Concatenate/combine all bitstring into a single bitstring
            BINCAT = strcat(BINSTR)
            // split the bitstring into groups of D matrix
            CHOP = strsplit(BINCAT)
            // Transform matrix CHOP to a 2D x D matrix
            BIN_ARR_STR = matrix(CHOP,[DIM 2*DIM])'
            // Convert the BIN_ARR_STR matrix values to intger
            B = strtod(BIN_ARR_STR)
            // convert 0 to -1 for easy use during eq. 2
            B(B==0) = -1
            // Eq.2
            S = repmat(P(point,:),[DIM*2 1]) + rand(DIM*2,DIM).*B.*(DELTA(point)/sqrt(2))
        end
        
        // Bounding variables that are out of bounds
        // Reinitialize VARIABLE VALUES that are out of bounds (less than 0; greater than 1)
        S(S<0) = rand(length(find(S<0)),1)
        S(S>1) = rand(length(find(S>1)),1)
        // scale to their actual variable values the 2p copies
        S = S.*(UPPER-LOWER)+LOWER
        SFIT = F1(S)   // evaluate
        
        square(LOWER,LOWER,UPPER,UPPER)
        plot(S(:,1),S(:,2),col(point))
        
        [SBESTFIT SBESTIND] = min(SFIT) // get best
        // return to scale the 2p copies
        S = (S-LOWER)./(UPPER-LOWER)
        // update Ppoint and it FITpoint with the fittest point from the generated copies
        P(point,:) = S(SBESTIND,:)  
        FIT(point) = SBESTFIT
    end
    
    // ARCHIVING
    if ARCHCOUNT==ARCHLIMIT // ARCH limit is reached
        if min(ARCHFIT)<FIT(1)
            // interchange P1 with the fittest point from the archive
            [ARCHBESTFIT ARCHIND] = min(ARCHFIT)
            P0 = ARCHBEST(ARCHIND,:)
            P(1,:) = ARCHBEST(ARCHIND,:)
            ARCHBEST(ARCHIND,:) = P0
            // do also for fitness
            FIT0 = ARCHBESTFIT
            FIT(1) = ARCHBESTFIT
            ARCHFIT(1) = FIT0
        end
        
        if min(ARCHFIT)<FIT(2)
            // make the fittest from the archive be the new P2
            [ARCHBESTFIT ARCHIND] = min(ARCHFIT)
            P(2,:) = ARCHBEST(ARCHIND,:)
            FIT(2) = ARCHBESTFIT
        end
        
        // update deltas/ radii
        DELTA(1) = DELTA(1)-(DELTA(1)/ALPHA)
        DELTA(2) = DELTA(2)+(DELTA(2)/ALPHA)
        
        // clear archive
        ARCHBEST = []
        ARCHFIT = []
        // set ARCHIVE COUNT to 0
        ARCHCOUNT = 0
        ARCHLIMIT = floor(rand()*(ARCHMAX+1-ARCHMIN))+ARCHMIN
    end
    
    // get ITERATION BEST
    [BESTFIT IND] = min(FIT)
    // update GBEST
    if BESTFIT<GBESTFIT
        GBESTFIT = BESTFIT
        GBESTPOS = P(IND,:).*(UPPER-LOWER)+LOWER
    end
    FITRUN(ITER) = GBESTFIT
    xtitle("ITER "+string(ITER))
    square(LOWER,LOWER,UPPER,UPPER)
    plot(P(:,1),P(:,2),'b.')
    plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
    xs2png(gcf(),'gif/yypo/ITER'+string(ITER)+'.png')
    clf()
end

//plot((1:MAXITER)',FITRUN,'g-')
