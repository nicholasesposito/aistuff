if (l_mandlvl and nlvinprof > 1):
    for j in range(1, nlv2wrt_tot + 1):
        jj = iord[j]
        nmNbtw = 0
        if (ibfms[drinfo_accum[1, jj]] != 0 and
            ibfms[drinfo_accum[2, jj]] != 0 and
            ibfms[drinfo_accum[3, jj]] != 0):
            nmNbtw = 1
            for k in range(j + 1, nlv2wrt_tot + 1):
                kk = iord[k]
                if (ibfms[drinfo_accum[1, kk]] != 0 and
                    ibfms[drinfo_accum[2, kk]] != 0 and
                    ibfms[drinfo_accum[3, kk]] != 0):
                    nmNbtw = nmNbtw + 1
                else:
                    break
            if (j <= 1):
                print()
                print('### PROBLEM - j <= 1 (=', j, ') in subr. ',
                      'sub2mem_mer, iord array underflow')
                print('              this should never happen')
                print()
                w3tage('PREPOBS_PREPACQC')
                errexit(61)
            jjm1 = iord[j - 1]
            jjpnmNbtw = iord[j + nmNbtw]
            pll = lvlsinprof[jjm1]
            pul = lvlsinprof[jjpnmNbtw]
            dtime_dlnp = (drinfo_accum[3, jjpnmNbtw] -
                          drinfo_accum[3, jjm1]) / log(pul / pll)
            lat_pul = drinfo_accum[2, jjpnmNbtw]
            lon_pul = drinfo_accum[1, jjpnmNbtw]
            lat_pll = drinfo_accum[2, jjm1]
            lon_pll = drinfo_accum[1, jjm1]
            if (int(lon_pul * 100) == int(lon_pll * 100)):
                dist_pul_pll = radius_e * abs(lat_pul - lat_pll) * deg2rad
            elif (int(lat_pul * 100) == int(lat_pll * 100)):
                dist_pul_pll = 2.0 * radius_e * asin(min(1.0, abs(cos(lat_pul * deg2rad) *
                                                                   sin((lon_pul - lon_pll) * 0.5 * deg2rad))))
            else:
                dist_pul_pll = 2.0 * radius_e * asin(min(1.0, sqrt(
                    (sin((lat_pul - lat_pll) * 0.5 * deg2rad)) ** 2
                    + cos(lat_pul * deg2rad) *
                    cos(lat_pll * deg2rad) *
                    (sin((lon_pul - lon_pll) * 0.5 * deg2rad)) ** 2
                )))
            if (int(drinfo_accum[3, jjpnmNbtw] * 100000) != int(drinfo_accum[3, jjm1] * 100000) and
                dist_pul_pll != 0):
                spd_pul_pll = dist_pul_pll / abs((drinfo_accum[3, jjpnmNbtw] -
                                                  drinfo_accum[3, jjm1]) * 3600)
                for k in range(nmNbtw):
                    jjpk = iord[j + k]
                    pml = lvlsinprof[jjpk]
                    drinfo_accum[3, jjpk] = drinfo_accum[3, jjm1] + dtime_dlnp * log(pml / pll)
                    dist2pml = spd_pul_pll * abs(drinfo_accum[3, jjpk] - drinfo_accum[3, jjm1]) * 3600
                    drinfo_accum[2, jjpk] = drinfo_accum[2, jjm1] + dist2pml / dist_pul_pll * (
                            drinfo_accum[2, jjpnmNbtw] - drinfo_accum[2, jjm1])
                    drinfo_accum[1, jjpk] = drinfo_accum[1, jjm1] + dist2pml / dist_pul_pll * (
                            drinfo_accum[1, jjpnmNbtw] - drinfo_accum[1, jjm1])
            else:
                delx = (drinfo_accum[1, jjpnmNbtw] -
                        drinfo_accum[1, jjm1]) / (nmNbtw + 1)
                dely = (drinfo_accum[2, jjpnmNbtw] -
                        drinfo_accum[2, jjm1]) / (nmNbtw + 1)
                for k in range(nmNbtw):
                    jjpk = iord[j + k]
                    pml = lvlsinprof[jjpk]
                    drinfo_accum[1, jjpk] = drinfo_accum[1, jjm1] + (k + 1) * delx
                    drinfo_accum[2, jjpk] = drinfo_accum[2, jjm1] + (k + 1) * dely
                    drinfo_accum[3, jjpk] = drinfo_accum[3, jjm1] + dtime_dlnp * log(pml / pll)
    jjmaxp = iord[nlv2wrt_tot]
    jjminp = iord[1]
    if (nlv2wrt_tot == 1):
        hdr2wrt[6] = 300 + (int(hdr2wrt[6]) % 100)
    elif (nlv2wrt_tot > 1 and
          (c_qc_accum[jjmaxp][10] == 'a' or
           c_qc_accum[jjmaxp][10] == 'A')):
        hdr2wrt[6] = 400 + (int(hdr2wrt[6]) % 100)
        hdr2wrt[2] = drinfo_accum[1, jjmaxp]
        hdr2wrt[3] = drinfo_accum[2, jjmaxp]
        hdr2wrt[4] = drinfo_accum[3, jjmaxp]
        hdr2wrt[5] = elv_accum[1, jjmaxp]
        hdr2wrt[12] = rpt_accum[1, jjmaxp]
        hdr2wrt[13] = tcor_accum[1, jjmaxp]
    elif (nlv2wrt_tot > 1 and
          (c_qc_accum[jjmaxp][10] == 'd' or
           c_qc_accum[jjmaxp][10] == 'D')):
        hdr2wrt[6] = 500 + (int(hdr2wrt[6]) % 100)
        hdr2wrt[2] = drinfo_accum[1, jjminp]
        hdr2wrt[3] = drinfo_accum[2, jjminp]
        hdr2wrt[4] = drinfo_accum[3, jjminp]
        hdr2wrt[5] = elv_accum[1, jjminp]
        hdr2wrt[12] = rpt_accum[1, jjminp]
        hdr2wrt[13] = tcor_accum[1, jjminp]
    hdr2wrt[10] = bmiss
    hdr2wrt[11] = bmiss



