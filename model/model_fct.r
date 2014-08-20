# model
#------------------------------------------
model = function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
	
	    #test
	       
#        TT = 0.2856791 
#        B=0.2057428 
#        M=0.3304669 
#        Ha=1590.129 
#        Hv=123.212

#        cat(TT, B, M, Ha, Hv, "\n")

		R=1-TT-B-M
#		Hv = Hv*dvv
#		Ha = Ha*dva
 #       cat("\n R; TT; B; M\n")
#        cat(R, TT, B, M)
		   
		# herbivores competition
		k = exp(-log(1/k0)*Hv/Ha)
		
		# moose and B
		mv_B = mvs/(1+(mvs/mv0-1)*(M+B))
		ma_B = mas/(1+(mas/ma0-1)*(M+B))
		
        # intakes
        if(Ha!=0) {
        Fa = ua*k*R/Ha
        Ga = k*(va*TT+wa*B+xa*M)/Ha
        }else{
        Fa=0
        Ga=0
        }
        if(Hv!=0) {
        Fv = uv*(1-k)*R/Hv
        Gv = (1-k)*(vv*TT+wv*B+xv*M)/Hv
        }else{
        Fv=0
        Gv=0
        } 
        
        Ia1 = taua*Fa/(mua+Fa)
        Ia2= (rhoa*Ga/(nua + Ga) )* (phia + (1-phia)/( 1 + exp(ra*( taua*Fa/(mua + Fa) -pa) ) ) )
        Iv1 = tauv*Fv/(muv+Fv)
        Iv2 = (rhov*Gv/(nuv + Gv) )* (phiv + (1-phiv)/( 1 + exp(rv*( tauv*Fv/(muv + Fv) -pv) ) ) )        
        # herbivore pressure
        
        if(R!=0) {
        PRB = (Ha*Ia1*kappaba)/(ua*R) + (Hv*Iv1*kappabv)/(uv*R)
        PRT = (Ha*Ia1*kappaba)/(ua*R) + (Hv*Iv1*kappabv)/(uv*R)
        }else{PRB = PRT =0}
        
        if(R!=1) 
        {
        sumOHa = omegatHa * T + omegabHa * B + (1-omegatHa-omegabHa) * M
        OtHa = omegatHa*T /sumOHa
        ObHa = omegabHa*B /sumOHa
        OmHa = (1-omegatHa-omegabHa)*M /sumOHa
 
        sumOHv = omegatHv * T + omegabHv * B + (1-omegatHv-omegabHv) * M
        OtHv = omegatHv*T /sumOHv
        ObHv = omegabHv*B /sumOHv
        OmHv = (1-omegatHv-omegabHv)*M /sumOHv
        
        if(TT!=0) {
        PTB = (OtHa*Ha*Ia2*kappaba)/(va*TT) + OtHv*Hv*Iv2*kappabv/(vv*TT)
        PTT = (OtHa*Ha*Ia2*kappata)/(va*TT) + OtHv*Hv*Iv2*kappatv/(vv*TT)
        }else{PTB = PTT =0}
        if(B!=0) {
        PBB = (ObHa*Ha*Ia2*kappaba)/(wa*B) + ObHv*Hv*Ia2*kappabv/(wv*B) 
        PBT = (ObHa*Ha*Ia2*kappata)/(wa*B) + ObHv*Hv*Ia2*kappatv/(wv*B) 
        }else{PBB = PBT =0}
        if(M!=0) {
        PMB = (OmHa*Ha*Ia2*kappaba)/(xa*M) + OmHv*Hv*Iv2*kappabv/(xv*M) 
        PMT = (OmHa*Ha*Ia2*kappata)/(xa*M) + OmHv*Hv*Iv2*kappatv/(xv*M) 
        }else{PMB = PMT =0}

        }else{
        PTB = PTT =0
        PBB = PBT =0
        PMB = PMT =0
        }
        
       
        # compute new alpha, theta and beta
        alphab_h = alphab*(1-PRB)
        alphat_h = alphat*(1-PRT)
        
        thetab_h = thetab*(1-PMB)
        thetat_h = thetat*(1-PMT)
        
        betab_h = (1/2)*(1 + betab * cos(pi*PTB) + (1-betab)*cos(pi*(1+PTT)) ) 
        betat_h = (1/2)*(1 + betat * cos(pi*PBT) + (1-betat)*cos(pi*(1+PBB)) )
        
	    # vegetation dynamic
	    R_T = alphat_h*(M+TT)*(1 - alphab_h*(M+B))
	    R_B = alphab_h*(M+B)*(1 - alphat_h*(M+TT))
	    R_M = alphat_h*(M+TT)*alphab_h*(M+B)
	    T_M = betab_h*(B+M)
	    B_M = betat_h*(TT+M)
        dTT = R*R_T + M*thetat_h - TT*(eps + T_M)
        dB = R*R_B + M*thetab_h - B*(eps + B_M)
        dM = R*R_M + TT*T_M + B*B_M - M*(eps + thetat_h + thetab_h)
                
        # herbivore dynamic
        dHa = Ha*( (1-ma_B)*((Ia1+Ia2)*ca-(pa+exp(-za*(Ia1 + Ia2)*ca))) - ma_B)
        dHv = Hv*( (1-mv_B) *((Iv1+Iv2)*cv-(pv+exp(-zv*(Iv1 + Iv2)*cv))) - mv_B)
        
        # negative values 
#        if(dHa<=(-Ha)){ dHa = round(-Ha)}
#        if(dHv<=(-Hv)) dHv = -Hv
#               

                      
		return(list(c(dTT, dB, dM, dHa, dHv)))
        
#        Ha_new = Ha+dHa
#        Hv_new = Hv+dHv
#        TT_new = TT+dTT
#        B_new = B + dB
#        M_new = M + dM
        
#        cat("\n dB; B; R; M; R_B; B_M\n")
#        cat(dB, B, R, M, R_B, B_M)

#        if(Ha_new<0){Ha_new=0}
#        if(Hv_new<0){Hv_new=0}
#        if(TT_new<0){TT_new=0}
#        if(B_new<0){B_new=0}
#        if(M_new<0){M_new=0}
        
#        list(c(TT_new, B_new, M_new, Ha_new, Hv_new))
	})
}

