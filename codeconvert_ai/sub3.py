rate_accum = bmiss
if nlvinprof == 0:
    print()
    print('### PROBLEM - into subr, sub2mem_mer with nlvinprof = ', 0)
    print('              this should never happen')
    print()
    call w3tage('PREPOBS_PREPACQC')
    call errexit(59)
orders(1, iwork, lvlsinprof, iord, nlvinprof, 1, lwr, 2)
nmandlvls = 0
nlv2wrt_tot = nlvinprof
if l_mandlvl and nlvinprof > 1:
    for i in range(1, maxmandlvls + 1):
        for j in range(1, nlvinprof + 1):
            jj = iord[j]
            jjp1 = iord[j + 1]
            if j < nlvinprof:
                if lvlsinprof[jj] < mandlvls[i] and lvlsinprof[jjp1] > mandlvls[i]:
                    if nlvinprof + nmandlvls + 1 > mxlv:
                        print(' #####> WARNING: THERE ARE MORE THAN ', mxlv, ' LEVELS IN THIS PROFILE -- WILL CONTINUE ON PROCESSING ONLY ', mxlv, ' LEVELS FOR THIS PROFILE')
                        cmxlv = str(mxlv)
                        call system('[ -n "$jlogfile" ] && $DATA/postmsg' + ' "$jlogfile" "***WARNING:' + cmxlv + ' AIRCRAFT ' + 'PROFILE LEVEL LIMIT EXCEEDED IN ' + 'PREPOBS_PREPACQC, ONLY ' + cmxlv + ' LEVELS ' + 'PROCESSED"')
                        break
                    nmandlvls = nmandlvls + 1
                    pll = lvlsinprof[jj]
                    pul = lvlsinprof[jjp1]
                    pqll = pevn_accum[2, jj, 1]
                    pqul = pevn_accum[2, jjp1, 1]
                    pml = mandlvls[i]
                    lvlsinprof[nlvinprof + nmandlvls] = mandlvls[i]
                    pevn_accum[1, nlvinprof + nmandlvls, 1] = pml / 10.
                    pevn_accum[2, nlvinprof + nmandlvls, 1] = max(pqll, pqul)
                    pevn_accum[3, nlvinprof + nmandlvls, 1] = nrlacqc_pc
                    pevn_accum[4, nlvinprof + nmandlvls, 1] = 98
                    cat_accum[1, nlvinprof + nmandlvls] = 7
                    if ibfms[tevn_accum[1, jj, 1]] == 0 and ibfms[tevn_accum[1, jjp1, 1]] == 0:
                        for iii in range(mxe4prof, 0, -1):
                            if ibfms[tevn_accum[1, jj, iii]] != 0:
                                nevents_t = iii
                            else:
                                nevents_t = iii
                                break
                        tll = tevn_accum[1, jj, nevents_t]
                        tqll = tevn_accum[2, jj, nevents_t]
                        for iii in range(mxe4prof, 0, -1):
                            if ibfms[tevn_accum[1, jjp1, iii]] != 0:
                                nevents_t = iii
                            else:
                                nevents_t = iii
                                break
                        tul = tevn_accum[1, jjp1, nevents_t]
                        tqul = tevn_accum[2, jjp1, nevents_t]
                        dt_dlnp = (tul - tll) / alog(pul / pll)
                        tml = tll + (dt_dlnp * (alog(pml / pll)))
                        tevn_accum[1, nlvinprof + nmandlvls, 1] = tml
                        tevn_accum[2, nlvinprof + nmandlvls, 1] = max(tqll, tqul)
                        tevn_accum[3, nlvinprof + nmandlvls, 1] = nrlacqc_pc
                        tevn_accum[4, nlvinprof + nmandlvls, 1] = 98
                    if ibfms[qevn_accum[1, jj, 1]] == 0 and ibfms[qevn_accum[1, jjp1, 1]] == 0:
                        for iii in range(mxe4prof, 0, -1):
                            if ibfms[qevn_accum[1, jj, iii]] != 0:
                                nevents_q = iii
                            else:
                                nevents_q = iii
                                break
                        qll = qevn_accum[1, jj, nevents_q]
                        qqll = qevn_accum[2, jj, nevents_q]
                        for iii in range(mxe4prof, 0, -1):
                            if ibfms[qevn_accum[1, jjp1, iii]] != 0:
                                nevents_q = iii
                            else:
                                nevents_q = iii
                                break
                        qul = qevn_accum[1, jjp1, nevents_q]
                        qqul = qevn_accum[2, jjp1, nevents_q]
                        dq_dlnp = (qul - qll) / alog(pul / pll)
                        qml = qll + (dq_dlnp * (alog(pml / pll)))
                        qevn_accum[1, nlvinprof + nmandlvls, 1] = qml
                        qevn_accum[2, nlvinprof + nmandlvls, 1] = max(qqll, qqul)
                        qevn_accum[3, nlvinprof + nmandlvls, 1] = nrlacqc_pc
                        qevn_accum[4, nlvinprof + nmandlvls, 1] = 98
                    else:
                        if ibfms[qbg_accum[2, jj]] == 0 and ibfms[qbg_accum[2, jjp1]] == 0:
                            qll = qbg_accum[2, jj]
                            qul = qbg_accum[2, jjp1]
                            dq_dlnp = (qul - qll) / alog(pul / pll)
                            qml = qll + (dq_dlnp * (alog(pml / pll)))
                            qbg_accum[2, nlvinprof + nmandlvls] = qml
                    if ibfms[zevn_accum[1, jj, 1]] == 0 and ibfms[zevn_accum[1, jjp1, 1]] == 0:
                        zll = zevn_accum[1, jj, 1]
                        zul = zevn_accum[1, jjp1, 1]
                        zqll = zevn_accum[2, jj, 1]
                        zqul = zevn_accum[2, jjp1, 1]
                        dz_dlnp = (zul - zll) / alog(pul / pll)
                        zml = zll + (dz_dlnp * (alog(pml / pll)))
                        zevn_accum[1, nlvinprof + nmandlvls, 1] = zml
                        zevn_accum[2, nlvinprof + nmandlvls, 1] = max(zqll, zqul)
                        zevn_accum[3, nlvinprof + nmandlvls, 1] = nrlacqc_pc
                        zevn_accum[4, nlvinprof + nmandlvls, 1] = 98
                    if ibfms[wuvevn_accum[1, jj, 1]] == 0 and ibfms[wuvevn_accum[1, jjp1, 1]] == 0 and ibfms[wuvevn_accum[2, jj, 1]] == 0 and ibfms[wuvevn_accum[2, jjp1, 1]] == 0:
                        for iii in range(mxe4prof, 0, -1):
                            if ibfms[wuvevn_accum[1, jj, iii]] != 0 or ibfms[wuvevn_accum[2, jj, iii]] != 0:
                                nevents_w = iii
                            else:
                                nevents_w = iii
                                break
                        ull = wuvevn_accum[1, jj, nevents_w]
                        vll = wuvevn_accum[2, jj, nevents_w]
                        uvqll = wuvevn_accum[3, jj, nevents_w]
                        for iii in range(mxe4prof, 0, -1):
                            if ibfms[wuvevn_accum[1, jjp1, iii]] != 0 or ibfms[wuvevn_accum[2, jjp1, iii]] != 0:
                                nevents_w = iii
                            else:
                                nevents_w = iii
                                break
                        uul = wuvevn_accum[1, jjp1, nevents_w]
                        vul = wuvevn_accum[2, jjp1, nevents_w]
                        uvqul = wuvevn_accum[3, jjp1, nevents_w]
                        du_dlnp = (uul - ull) / alog(pul / pll)
                        dv_dlnp = (vul - vll) / alog(pul / pll)
                        uml = ull + (du_dlnp * (alog(pml / pll)))
                        vml = vll + (dv_dlnp * (alog(pml / pll)))
                        wuvevn_accum[1, nlvinprof + nmandlvls, 1] = uml
                        wuvevn_accum[2, nlvinprof + nmandlvls, 1] = vml
                        wuvevn_accum[3, nlvinprof + nmandlvls, 1] = max(uvqll, uvqul)
                        wuvevn_accum[4, nlvinprof + nmandlvls, 1] = nrlacqc_pc
                        wuvevn_accum[5, nlvinprof + nmandlvls, 1] = 98
    nlv2wrt_tot = nlvinprof + nmandlvls
    orders(1, iwork, lvlsinprof, iord, nlv2wrt_tot, 1, lwr, 2)
write(41, *) 'nlv2wrt_tot=', nlv2wrt_tot, 'c_acftreg=', c_acftreg1
err_tspline = 0



