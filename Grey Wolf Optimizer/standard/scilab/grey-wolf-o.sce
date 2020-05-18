
clear
clc

exec('objective.sce');
disp('RUNNING...')
rand('seed',getdate('s'))

POPSIZE = 50
DIM = 2
UPPER = 100
LOWER = -100
MAXITER = 50

// INITIALIZE WOLF POPULATION
WOLF = rand(POPSIZE,DIM).*(UPPER-LOWER)+LOWER
// EVALUATE WOLVES FITNESS
FIT = F1(WOLF)
// SORT WOLF ACCORDING TO FITNESS
[FIT SORTIND] = gsort(FIT,'lr','i')
WOLF = WOLF(SORTIND,:)
// GET GBEST
GBESTPOS = WOLF(1,:)
GBESTFIT = FIT(1)
xtitle("INITIALIZATION")
square(LOWER,LOWER,UPPER,UPPER)
plot(WOLF(:,1),WOLF(:,2),'b.')
plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
xs2png(gcf(),'gif/gwo/ITER0.png')
clf()

for ITER=1:MAXITER
    a = 2 - (2-0)*((ITER-1)/MAXITER)
    A1 = rand(POPSIZE,DIM).*(2*a) - a  // Eq. 3.3
    A2 = rand(POPSIZE,DIM).*(2*a) - a  // Eq. 3.3
    A3 = rand(POPSIZE,DIM).*(2*a) - a  // Eq. 3.3
    C1 = rand(POPSIZE,DIM).*a  // Eq. 3.4
    C2 = rand(POPSIZE,DIM).*a  // Eq. 3.4
    C3 = rand(POPSIZE,DIM).*a  // Eq. 3.4
    D_ALPHA = abs(repmat(WOLF(1,:),[POPSIZE 1]).*C1 - WOLF) // Eq.3.5 part 1
    D_BETA = abs(repmat(WOLF(2,:),[POPSIZE 1]).*C2 - WOLF) // Eq.3.5 part 1
    D_DELTA = abs(repmat(WOLF(3,:),[POPSIZE 1]).*C3 - WOLF) // Eq.3.5 part 1
    X1 = repmat(WOLF(1,:),[POPSIZE 1]) - A1.*D_ALPHA // Eq.3.6 part 2
    X2 = repmat(WOLF(2,:),[POPSIZE 1]) - A2.*D_BETA // Eq.3.6 part 2
    X3 = repmat(WOLF(3,:),[POPSIZE 1]) - A3.*D_DELTA// Eq.3.6 part 2
    WOLF = (X1 + X2 + X3)./3
    // AMMEND OUT OF BOUNDS
    WOLF = ammend(WOLF,UPPER,LOWER)
    //WOLF = penalty(WOLF,ubx,lbx,lby,uby)
    // EVALUATE WOLVES FITNESS
    FIT = F1(WOLF)
    // SORT WOLF ACCORDING TO FITNESS
    [FIT SORTIND] = gsort(FIT,'lr','i')
    WOLF = WOLF(SORTIND,:)
    // UPDATE GBEST
    if FIT(1) < GBESTFIT
        GBESTFIT = FIT(1)
        GBESTPOS = WOLF(1,:)
    end
    FITRUN(ITER) = GBESTFIT
    
    xtitle("ITER "+string(ITER))
    square(LOWER,LOWER,UPPER,UPPER)
    plot(WOLF(:,1),WOLF(:,2),'b.')
    plot(GBESTPOS(:,1),GBESTPOS(:,2),'g.')
    xs2png(gcf(),'gif/gwo/ITER'+string(ITER)+'.png')
    clf()
    
end
plot((1:MAXITER)',GBEST,'g-')
