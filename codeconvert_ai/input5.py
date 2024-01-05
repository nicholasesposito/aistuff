def slalom_uspline(n, bend, hgxn, yn, off, halfgate, q, ya, en, ita, ittot, FF):
    nita = 50
    nitb = 80
    xs = hgxn * halfgate
    hgxp = hgxn[0] - 1
    on = [False] * n
    FF = False
    for i in range(n):
        if off[i]:
            on[i] = False
            continue
        on[i] = hgxn[i] > hgxp
        if on[i]:
            hgxp = hgxn[i]
    ittot = 1
    fit_guspline(n, xs, yn, on, qt, jump, yat, en, FF)
    ena = en
    if FF:
        print('In slalom_uspline; failure flag raised in call to fit_guspline')
        print('at initialization of A loop')
        return
    for ita in range(1, nita + 1):
        q = qt
        ya = yat
        j = 0
        k = 0
        sjmin = 0
        for i in range(n):
            if not on[i]:
                continue
            sj = -bend[i] * jump[i]
            if sj < 0:
                j = i
                on[i] = False
            else:
                k += 1
        if j == 0:
            break
        if k == 0:
            on[j] = True
        for itb in range(1, nitb + 1):
            fit_guspline(n, xs, yn, on, qt, jump, yat, en, FF)
            if FF:
                print('In slalom_uspline; failure flag raised in call to fit_guspline')
                print('at B loop, iterations ita,itb = ', ita, itb)
                return
            ittot += 1
            fit_uspline()
            j = 0
            sjmin = u1
            for i in range(n):
                if on[i] or off[i]:
                    continue
                sj = bend[i] * (yn[i] - yat[i])
                if sj < 0:
                    sj = (yn[i] - ya[i]) / (yat[i] - ya[i])
                    if sj < sjmin:
                        j = i
                        sjmin = sj
            if j == 0:
                break
            solution as A
            ya = ya + sjmin * (yat - ya)
            q = q + sjmin * (qt - q)
            on[j] = True
        maxitb = max(maxitb, itb)
        if itb > nitb:
            FF = True
            print('In slalom_uspline; exceeding the allocation of B iterations')
            return
        q = qt
        ya = yat
        if en >= ena:
            print('In slalom_uspline; energy failed to decrease')
            break
        ena = en
    if ita > nita:
        FF = True
        print('In slalom_uspline; exceeding the allocation of A iterations')
        return

def convertd(n, halfgate, tdata, hdata, phof, doru, idx, hgts, hs, descending, FF):
    hour = 3600
    FF = False
    if len(hdata) != n:
        raise ValueError('In convertd; inconsistent dimensions of hdata')
    if len(tdata) != n:
        raise ValueError('In convertd; inconsistent dimensions of tdata')
    if len(hs) != n:
        raise ValueError('In convertd; inconsistent dimensions of hs')
    if len(hgts) != n:
        raise ValueError('In convertd; inconsistent dimensions of hgts')
    hs = hdata
    upsign = 0
    gate = halfgate * 2
    for i in range(n):
        hgts[i] = 2 * round(tdata[i] * hour / gate)
        if phof[i] == 5:
            upsign = 1
        if phof[i] == 6:
            upsign = -1
    doru = 0
    if upsign > 0:
        doru = 2
    else:
        doru = 1
    if n == 1:
        return
    if hgts[0] > hgts[n - 1]:
        for i in range(n // 2):
            j = n - i - 1
            hgs = hgts[i]
            hgts[i] = hgts[j]
            hgts[j] = hgs
            s = hs[i]
            hs[i] = hs[j]
            hs[j] = s
    if upsign == 1:
        descending = False
    elif upsign == -1:
        descending = True
    else:
        descending = hs[n - 1] < hs[0]
        if descending:
            upsign = -1
            print('mainly DESCENDING')
        else:
            upsign = 1
            print('mainly ASCENDING')
    idx = list(range(1, n + 1))
    for i in range(1, n):
        for ii in range(i):
            if hgts[i] < hgts[ii] or (hgts[i] == hgts[ii] and upsign * (hs[i] - hs[ii])) < u0:
                hgs = hgts[i]
                hgts[i] = hgts[ii]
                hgts[ii] = hgs
                s = hs[i]
                hs[i] = hs[ii]
                hs[ii] = s
                hgts
                j = idx[i]
                idx[i] = idx[ii]
                idx[ii] = j
    for i in range(1, n):
        if hgts[i] < hgts[i - 1]:
            print('In convertd; time sequence not monotonic', i, hgts[i], hgts[i - 1])
            FF = True
            return
    for i in range(1, n):
        if upsign * (hs[i] - hs[i - 1]) < u0:
            print('In convertd; height sequence not monotonic')
            FF = True
            return

def convertd_back(n, halfgate, wdata, tdata, ws, hgts, idx, descending):
    if len(wdata) != n:
        raise ValueError('In convertd; inconsistent dimensions of wdata')
    if len(tdata) != n:
        raise ValueError('In convertd; inconsistent dimensions of tdata')
    if len(ws) != n:
        raise ValueError('In convertd; inconsistent dimensions of ws')
    if len(hgts) != n:
        raise ValueError('In convertd; inconsistent dimensions of hgts')
    wn = [0] * n
    hgtn = [0] * n
    for i in range(n):
        ii = idx[i]
        hgtn[ii] = hgts[i]
        wn[ii] = ws[i]
    wdata = wn
    tdata = hgtn * halfgate
    if descending or n == 1:
        return
    for i in range(n // 2):
        j = n - i - 1
        s = tdata[i]
        tdata[i] = tdata[j]
        tdata[j] = s
        s = wdata[i]
        wdata[i] = wdata[j]
        wdata[j] = s



