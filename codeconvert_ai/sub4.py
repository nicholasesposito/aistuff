if ((nlv2wrt_tot > 1) and tsplines):
    nh = 0
    for j in range(1, nlv2wrt_tot + 1):
        jj = iord[j]
        if (ibfms[drinfo_accum[3, jj]] == 0):
            nh += 1
    nh2 = nh * 2
    halfgate = 30.0
    print('halfgate=', halfgate)
    idx = [0] * nh
    pof = [0] * nh
    tdata = [0] * nh
    hdata = [0] * nh
    wdata = [0] * nh
    te = [0] * nh
    hgts = [0] * nh
    hs = [0] * nh
    dhdt = [0] * nh
    maxita = 0
    maxitb = 0
    maxrts = 0
    maxit = 0
    nh = 0
    for j in range(1, nlv2wrt_tot + 1):
        jj = iord[j]
        if (ibfms[drinfo_accum[3, jj]] == 0):
            nh += 1
            tdata[nh] = drinfo_accum[3, jj]
            hdata[nh] = zevn_accum[1, jj, 1]
            pof[nh] = int(acft_seq_accum[2, jj])
            print('tdata,hdata,pof=', nh, tdata[nh], hdata[nh], pof[nh])
    convertd(nh, halfgate, tdata, hdata, pof, doru, idx, hgts, hs, descending, FF)
    if (FF):
        print("WARNING: tspline err in utility pspl, coming out of subr. convertd - use finite difference method")
        print("WARNING: tspline err in utility pspl, coming out of subr. convertd - use finite difference method")
        err_tspline = 1
        goto 666
    if (descending):
        print('set descending')
    else:
        print('set ascending')
    count_gates(nh, hgts[0:nh], mh)
    m = mh * 2
    hgtp = [0] * m
    hp = [0] * m
    qbest = [0] * m
    habest = [0] * m
    modebest = [0] * mh
    best_slalom(nh, mh, doru, hgts, hs, halfgate, bigT, hgtp, hp, qbest, habest, enbest, modebest, maxita, maxitb, maxit, maxrts, FF)
    print('maxita,maxitb,maxit,maxrts=', maxita, maxitb, maxit, maxrts)
    if (FF):
        print("WARNING: tspline err in utility pspl, coming out of subr. best_slalom - use finite difference method")
        print("WARNING: tspline err in utility pspl, coming out of subr. best_slalom - use finite difference method")
        err_tspline = 1
        goto 666
    bnewton(nh, m, bigT, halfgate, hgts, hs, hgtp, habest, qbest, te[0:nh], dhdt[0:nh], FF)
    if (FF):
        print("WARNING: tspline err in utility pspl, coming out of subr. bnewton - use finite difference method")
        print("WARNING: tspline err in utility pspl, coming out of subr. bnewton - use finite difference method")
        err_tspline = 1
        goto 666
    convertd_back(nh, halfgate, wdata, tdata, dhdt, hgts, idx, descending)
    for j in range(1, nh + 1):
        print('hgts,hs,dhdt,wdata=', j, hgts[j], hs[j], dhdt[j], wdata[j])
    nh = 0
    for j in range(1, nlv2wrt_tot + 1):
        jj = iord[j]
        if (ibfms[drinfo_accum[3, jj]] == 0):
            nh += 1
            rate_accum[jj] = wdata[nh]
            print('j,z,rate=', j, zevn_accum[1, jj, 1], rate_accum[jj])
    continue_666
    if (allocated(idx)):
        del idx
    if (allocated(pof)):
        del pof
    if (allocated(tdata)):
        del tdata
    if (allocated(hdata)):
        del hdata
    if (allocated(wdata)):
        del wdata
    if (allocated(te)):
        del te
    if (allocated(hgts)):
        del hgts
    if (allocated(hs)):
        del hs
    if (allocated(dhdt)):
        del dhdt
    if (allocated(hgtp)):
        del hgtp
    if (allocated(hp)):
        del hp
    if (allocated(qbest)):
        del qbest
    if (allocated(habest)):
        del habest
    if (allocated(modebest)):
        del modebest
endif
if (((nlv2wrt_tot > 1) and (not tsplines)) or (err_tspline > 0)):
    for j in range(1, nlv2wrt_tot + 1):
        jj = iord[j]
        print('j,ord,z,t,pof=', j, jj, zevn_accum[1, jj, 1], drinfo_accum[3, jj], acft_seq_accum[1, jj], acft_seq_accum[2, jj])
    for j in range(1, nlv2wrt_tot + 1):
        jj = iord[j]
        jkp = 0
        jkm = 0
        jjp1 = 0
        jjm1 = 0
        if (j == nlv2wrt_tot):
            if (ibfms[drinfo_accum[3, jj]] == 0):
                jjp1 = jj
                jkp = j
        else:
            for jk in range(j + 1, nlv2wrt_tot + 1):
                jjp = iord[jk]
                if (jjp > nlvinprof):
                    continue
                if (ibfms[drinfo_accum[3, jjp]] == 0):
                    jjp1 = jjp
                    jkp = jk
                    break
        if (j == 1):
            if (ibfms[drinfo_accum[3, jj]] == 0):
                jjm1 = jj
                jkm = j
        else:
            for jk in range(j - 1, 0, -1):
                jjm = iord[jk]
                if (jjm > nlvinprof):
                    continue
                if (ibfms[drinfo_accum[3, jjm]] == 0):
                    jjm1 = jjm
                    jkm = jk
                    break
        if ((jjp1 != 0) and (jjm1 != 0)):
            dt = (drinfo_accum[3, jjp1] - drinfo_accum[3, jjm1]) * 3600.
            c1_jk = 0
            c2_jk = 0
            while ((abs(dt) < 60.) and ((jkp + c1_jk <= nlv2wrt_tot) or (jkm - c2_jk >= 1))):
                jjp2 = 0
                jjm2 = 0
                c1_jk += 1
                c2_jk += 1
                dt_new = dt
                while (jkp + c1_jk <= nlv2wrt_tot and iord[jkp + c1_jk] > nlvinprof):
                    c1_jk += 1
                if (jkp + c1_jk <= nlv2wrt_tot and iord[jkp + c1_jk] <= nlvinprof):
                    jjp = iord[jkp + c1_jk]
                    if (ibfms[drinfo_accum[3, jjp]] == 0):
                        jjp2 = jjp
                        dt_new = (drinfo_accum[3, jjp2] - drinfo_accum[3, jjm1]) * 3600.
                if (abs(dt_new) >= 60.):
                    if (jjp2 != 0):
                        jjp1 = jjp2
                    break
                while (jkm - c2_jk >= 1 and iord[jkm - c2_jk] > nlvinprof):
                    c2_jk += 1
                if (jkm - c2_jk >= 1 and iord[jkm - c2_jk] <= nlvinprof):
                    jjm = iord[jkm - c2_jk]
                    if (ibfms[drinfo_accum[3, jjm]] == 0):
                        jjm2 = jjm
                        dt_new = (drinfo_accum[3, jjp1] - drinfo_accum[3, jjm2]) * 3600.
                if (abs(dt_new) >= 60.):
                    if (jjm2 != 0):
                        jjm1 = jjm2
                    break
                if ((jjp2 != 0) and (jjm2 != 0)):
                    dt_new = (drinfo_accum[3, jjp2] - drinfo_accum[3, jjm2]) * 3600.
                    if (abs(dt_new) >= 60.):
                        if (jjp2 != 0):
                            jjp1 = jjp2
                        if (jjm2 != 0):
                            jjm1 = jjm2
                        break
            dt = (drinfo_accum[3, jjp1] - drinfo_accum[3, jjm1]) * 3600.
            zul = zevn_accum[1, jjp1, 1]
            zll = zevn_accum[1, jjm1, 1]
            if (abs(dt) > 0.):
                rate_accum[jj] = (zul - zll) / dt
            print('fj,dt,rate_accum=', j, dt, rate_accum[jj])
            print('')



