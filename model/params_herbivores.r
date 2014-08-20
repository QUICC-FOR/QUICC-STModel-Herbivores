#    size = 200000 #ha
    u = 0.5*1000; v=200; w = 100; x=200
	#params fixes
	k0 = 0.5

    omegatHa = 0.3
    omegabHa = 0.35

    omegatHv = 0.325
    omegabHv = 0.35

    kappata = 0.25
    kappaba = 0.07

    kappatv = 0.265
    kappabv = 0.175
    
    # moose
    dva = 40000
    ua = u*dva ; va = u*dva ;  wa= u*dva ; xa= u*dva

    za = 1
    ca = 0.79
    ma0 = 0.25
    mas = 0.25
    taua = 13
    mua = (taua/0.07)*(dva/100)/sqrt(314)
    rhoa = 6
    nua = (rhoa/0.15)*(dva/100)/sqrt(314)
    phia = 0.5
    pa = 1.7
    ra = 5

    #deer
    dvv = 1000 # ha
    uv = u*dvv ; vv = u*dvv ;  wv= u*dvv ; xv= u*dvv


    zv = 1
    cv = 0.5
    mv0 = 0.2
    mvs = 0.5
    tauv = 17.5
    muv = (tauv/0.07)*(dvv/100)/sqrt(518)
    rhov = 4.5
    nuv = (rhov/0.15)*(dvv/100)/sqrt(518)
    pv = 1.87
    phiv = 0.57
    rv = 5

