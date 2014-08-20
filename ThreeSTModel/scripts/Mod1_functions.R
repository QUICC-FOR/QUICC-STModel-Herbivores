rm(list=ls())
# Range kutta 
#------------------------------------------
library(deSolve)

#------------------------------------------

#fixed params
alpha0=0.07 # succession rate (1 / age maturity)
delta=0.05 # disturbance rate (intensity * frequency = 0.75 * 1/25)
f=1 # trees fecundity

taug = 6; thetag = 0
mug = 20; nug = 15
rg=5; pg = mg = 0.2;

taub = 6; thetab=5
mub = 18; nub = 25
rb=5; pb = mb = 0.2;
k0=0.5

#______________________________________________________________________

run <- function(Hb, Hg, comp , k0=0.5, u=1000, w=4000)
{

# initial conditions
T0 = c(T=0.5, S=0.5, Hb=Hb, Hg = Hg)

# simulation conditions
# parameters

comp0= comp # competitiveness of seedlings on grasses when herbivore is absent




#------------------------------------------
model = function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
	
	    
	    G=1-T-S        

        if(Hg==0) k=1;
        if(Hb==0) k=0;
        if(Hg*Hg!=0) k = exp(-log(1/k0)*Hg/Hb)

        # intakes
        if(Hb!=0) {
        Fb = u*k*S/Hb
        Gb = k*w*G/Hb
        }else{
        Fb=0
        Gb=0
        }
        if(Hg!=0) {
        Fg = (1-k)*w*G/Hg
        Gg = u*(1-k)*S/Hg
        }else{
        Fg=0
        Gg=0
        } 
        
        Ib1 = taub*Fb/(mub+Fb)
        Ib2 = (thetab*Gb/(nub + Gb) )* (1/( 1 + exp(rb*( taub*Fb/(mub + Fb) -pb) ) ) )
        Ig1 = taug*Fg/(mug+Fg)
        Ig2 = (thetag*Gg/(nug + Gg) )* (1/( 1 + exp(rg*( taug*Fg/(mug + Fg) -pg) ) ) )        
        # herbivore pressure
  
        if(S!=0) {
        PS = (Hb*Ib1 + Hg*Ig2)/(u*S)
        }else{PS =0}
        if(G!=0) {
        PG = (Hb*Ib2 + Hb*Ig1)/(w*G)
        }else{PG =0}
        
        # compute alpha and c
        alpha = alpha0*(1-PS)
        comp = PG - comp0*(PG + PS) + comp0
        
	    # vegetation model
	    dT= alpha*S - delta*T
		dS=f*T*comp*G - alpha*S
		
		#herbivore model
		Ib = Ib1 + Ib2
		Ig = Ig1 + Ig2
		dHb = Hb*(Ib-(pb+exp(-Ib)+mb))
		dHg = Hg*(Ig-(pg+exp(-Ig)+mg))
	
		return(list(c(dT, dS, dHb, dHg)))
	})
}

pars = c(alpha0=alpha0, 
delta=delta,
f=f, 
comp0=comp0,
k0=k0,
u=u,
w=w,
taug=taug,
thetag=thetag,
mug=mug,
nug=nug,
rg=rg,
pg=pg,
mg=mg,
taub=taub,
thetab=thetab,
mub=mub,
nub=nub,
rb=rb,
pb=pb,
mb=mb
)

#------------------------------------------
# run simulation
out = stode(func=model, y=T0, parms = pars, positive = TRUE)
#summary(out)

eq = out$y
eq["G"] = 1-eq["T"] -eq["S"]
return(eq[c(1:2,5, 3, 4)])
}

#------------------------------------------

runNoH <- function(env, u0=1000, w0=4000)
{

## environmental variation
comp0 = env[1]
u = env[2]*u0
w = env[2]*w0

eq = run(Hb= 0, Hg = 0, comp = comp0, u = u , w = w)

return(eq)
}

#------------------------------------------

runHb <- function(env, u0=1000, w0=4000)
{

## environmental variation
comp0 = env[1]
u = env[2]*u0
w = env[2]*w0

eq = run(Hb= 1000, Hg = 0, comp = comp0, u = u , w = w)

return(eq)
}

#-----------------------------------------------------------
#-----------------------------------------------------------


runHg <- function(env, u0=1000, w0=4000)
{
## environmental variation
comp0 = env[1]
u = env[2]*u0
w = env[2]*w0

eq = run(Hb= 0, Hg = 1000, comp = comp0, u = u , w = w)


return(eq)
}

#-----------------------------------------------------------
#-----------------------------------------------------------


runHH <- function(env, u0=1000, w0=4000)
{

comp0 = env[1]
u = env[2]*u0
w = env[2]*w0

eq = run(Hb= 1000, Hg = 1000, comp = comp0, u = u , w = w)


return(eq)
}


##
#
##------------------------------------------
## coexistence
##------------------------------------------
library(pracma) #lambertWp
library(rootSolve)

coexistence <- function(GTS, u=1000, w = 4000)
{


# carrying capacities 
#--------

Ikb = mb + pb + lambertWp(exp(-(mb+pb)))
Ikg = mg + pg + lambertWp(exp(-(mg+pg)))


# H equilibrium functions
#----------

# to solve Hb for Ib - Ikb = 0
Ib_Ikb = function(Hb, S, G, Hg)
{
if(Hb!=0){
    if(Hg==0){
    k = 1
    }else k = exp(-log(1/k0)*Hg/Hb)
    
    Ab = taub *k*u*S/(Hb*mub +k*u*S)
    res = Ab + ( thetab *w*G/(Hb*nub +w*G) )* (1/ (1+exp( rb*( Ab -pb-mb ))) ) -Ikb  
}else res=0;

return(res)
}


# to solve Hb for Ig - Ikg =0
Ig_Ikg = function(Hb, S, G, Hg)
{
if(Hg!=0){
    if(Hb==0){
    k = 0
    }else k = exp(-log(1/k0)*Hg/Hb)
    
    Ag = taug *w*G/(Hg*mug +w*G)
    res = Ag + ( thetag *(1-k)*u*S/(Hg*nug +(1-k)*u*S) )* (1/ (1+exp( rg*( Ag -pg-mg ))) ) -Ikg  
}else res=0;

return(res)
}


# to solve Hg for Ig - Ikg =0
Ig_Ikg2 = function(Hg, S, G, Hb)
{
if(Hg!=0){
    if(Hb==0){
    k = 0
    }else k = exp(-log(1/k0)*Hg/Hb)
    
    Ag = taug *w*G/(Hg*mug +w*G)
    res = Ag + ( thetag *(1-k)*u*S/(Hg*nug +(1-k)*u*S) )* (1/ (1+exp( rg*( Ag -pg-mg ))) ) -Ikg  
}else res=0;

return(res)
}



eq_fct = function(Hg, GTS)
{

#GTS = c(.5, .5, 0) ; Hg=0; # debug
G = as.numeric(GTS[1])
T = as.numeric(GTS[2])
S = as.numeric(GTS[3])

# solve Ib - Ikb = 0 pour Hb
# cherche solution positive sinon 0 est une solution quelque-soit Hg
hbmin = as.numeric(Ib_Ikb(0.0001, S=S, G=G, Hg=Hg))
hbmax = as.numeric(Ib_Ikb(100000, S=S, G=G, Hg=Hg))
if(hbmin*hbmax<0) 
{
Hb_eq1 = uniroot(Ib_Ikb, S=S, G=G, Hg=Hg, interval=c(0.0001, 100000), f.lower = hbmin, f.upper = hbmax, tol = 0.001)$root
} else Hb_eq1 = 0;


# solve Ig -Ikg = 0 pour Hb
hbmin = as.numeric(Ig_Ikg(0.0001, S=S, G=G, Hg=Hg))
hbmax = as.numeric(Ig_Ikg(100000, S=S, G=G, Hg=Hg))
if(hbmin*hbmax<0) 
{
Hb_eq2 = uniroot(Ig_Ikg, S=S, G=G, Hg=Hg, interval=c(0.0001, 100000), f.lower = hbmin, f.upper = hbmax, tol = 0.001)$root
} else Hb_eq2 = Hb_eq1;

# res should be 0 to solve both Ikb-Ib = 0 and Ikg-Ig = 0
res = Hb_eq1 - Hb_eq2

return(res)
}

### solve equilibrium function
eq_H = function(GTS)
{

#GTS = c(0.9, 0, .1) #debug

G= as.numeric(GTS[1])
S=as.numeric(GTS[3])
T=as.numeric(GTS[2])


# need to find the minimun to have an interval around zero
Hgx = seq(0.0001, 100000, 10)
res = unlist(lapply(Hgx, eq_fct, GTS=as.numeric(GTS[1:3])))
#plot(Hgx, res, type="b")
#abline(h=0)

# special case (we will check them later if necessary)
if(length(na.omit(res))==1)res=NA

# test if there is possible equilibrium of both
if(sum(is.na(res))!=length(res))
{
print(GTS)
if(max(res, na.rm=TRUE)>0 & min(res, na.rm=TRUE)<0 )
{

    imin = Hgx[which.min(res)]
    imax = Hgx[which.max(res)]
    hgmin = eq_fct(imin, as.numeric(GTS[1:3]))
    hgmax = eq_fct(imax, as.numeric(GTS[1:3]))
    if(hgmin*hgmax<0) 
    {
    Hg_eq = uniroot(eq_fct, GTS=as.numeric(GTS[1:3]), interval=c(imin, imax))$root
    } else {
    if(eq_fct(0, GTS=as.numeric(GTS[1:3]))==0) {
    Hg_eq = 0
    }else Hg_eq=NA
    }

   # find Hb back pour Hg_eq
    if(!is.na(Hg_eq))
    {

    hbmin = as.numeric(Ib_Ikb(0.0001, S=S, G=G, Hg=Hg_eq))
    hbmax = as.numeric(Ib_Ikb(100000, S=S, G=G, Hg=Hg_eq))
    if(hbmin*hbmax<0) 
    {
    Hb_eq = uniroot(Ib_Ikb, S=S, G=G, Hg=Hg_eq, interval=c(0.0001, 100000), f.lower = hbmin, f.upper = hbmax, tol = 0.001)$root
    } else Hb_eq = 0;

    }else {
    print("pb:")
    print(GTS)
    Hb_eq = NA
    }

}else{
## res == 0 pour tout Hg>0, which mean Hb=0
Hb_eq= 0

   # find Hg back pour Hb=0
    hgmin = as.numeric(Ig_Ikg2(0.0001, S=S, G=G, Hb=0))
    hgmax = as.numeric(Ig_Ikg2(100000, S=S, G=G, Hb=0))
    if(hgmin*hgmax<0) 
    {
    Hg_eq = uniroot(Ig_Ikg2, S=S, G=G, Hb=0, interval=c(0.0001, 100000), f.lower = hgmin, f.upper = hgmax, tol = 0.001)$root
    } else Hg_eq = 0;



}
}else {
# res = NA pour tout Hg (eventuellement sauf 1): NE devrait pas arriver
Hb_eq= Hg_eq = NA;
}
return(c(Hb_eq, Hg_eq))
}

##------------------
return(eq_H(GTS))

}
