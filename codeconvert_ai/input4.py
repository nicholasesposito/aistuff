def set_gates(nh, mh, doru, hgts, hs, hgtn, hn, code, FF):
    FF = False
    n = nh * 2
    hgtp = hgts[0] - 1
    imh = 0
    for i in range(nh):
        i2 = i * 2
        i2m = i2 - 1
        hp = hs[i]
        if hgts[i] > hgtp:
            imh += 1
            hgtp = hgts[i]
            hgtn[0][imh-1] = hgtp - 1
            hgtn[1][imh-1] = hgtp + 1
            hn[:,:,imh-1] = hp
        elif hgts[i] < hgtp:
            FF = True
            print('In set_gates; data are not temporally monotonic')
            return
        else:
            hn[0][0][imh-1] = max(hn[0][0][imh-1], hp)
            hn[0][1][imh-1] = min(hn[0][1][imh-1], hp)
            hn[1][1][imh-1] = max(hn[1][1][imh-1], hp)
            hn[1][0][imh-1] = min(hn[1][0][imh-1], hp)
    if imh != mh:
        raise ValueError('In set_gates; inconsistent gate tallies, imh and mh')
    if mh == 1:
        code[0] = 4 * doru
        return

def set_posts(mh, mode, hgtn, hn, bend, hgtp, hp, off):
    off = False
    for i in range(mh):
        im = i - 1
        modei = mode[i]
        i2 = i * 2
        i2m = i2 - 1
        i2mm = i2 - 2
        hgtp[i2m] = hgtn[0][i]
        hgtp[i2] = hgtn[1][i]
        hp[i2m] = hn[0][modei][i]
        hp[i2] = hn[1][modei][i]
        if i > 1:
            if hgtprev == hgtp[i2m]:
                if hprev == hp[i2m]:
                    off[i2m] = True
                if mode[im] == 2 and modei == 1:
                    if hprev <= hp[i2m]:
                        off[i2mm] = True
                    else:
                        off[i2m] = True
                elif mode[im] == 1 and modei == 2:
                    if hprev <= hp[i2m]:
                        off[i2m] = True
                    else:
                        off[i2mm] = True
        bend[i2m] = modei * 2 - 3
        bend[i2] = -bend[i2m]
        hgtprev = hgtp[i2]
        hprev = hp[i2]

def count_routes(n, code, count, FF):
    FF = False
    flag = True
    count = 0
    while True:
        next_route(n, code, mode, flag)
        if flag:
            return
        count += 1
    FF = (count > ihu)
    if FF:
        print('In count_routes; number of routes exceeds allowance =', ihu)

def list_routes(n, code):
    print('List all route combinations of', n, 'allowed passage modes')
    flag = True
    for i in range(ihu):
        next_route(n, code, mode, flag)
        if flag:
            print('In list_routes; List of routes complete')
            break
        print(i, mode)

def next_route(n, code, mode, flag):
    options = [[0, 1, 2], [0, 1, 2], [0, 1, 2], [0, 0, 0], [1, 1, 1], [2, 2, 2]]
    firstmode = [1, 1, 2]
    modeim = 1
    if flag:
        for i in range(n):
            option = options[code[i]][modeim]
            mode[i] = firstmode[option]
            modeim = mode[i]
        flag = False
        return
    for i in range(n, 0, -1):
        im = i - 1
        if i > 1:
            modeim = mode[im]
        else:
            modeim = 1
        option = options[code[i-1]][modeim]
        if option > 0 or mode[i-1] == 2:
            continue
        mode[i-1] = 2
        modejm = mode[i-1]
        for j in range(i+1, n+1):
            option = options[code[j-1]][modejm]
            mode[j-1] = firstmode[option]
            modejm = mode[j-1]
        return
    flag = True

def slalom_tspline(n, bend, hgxn, yn, off, bigX, q, ya, en, ita, maxitb, ittot, FF):
    xs = hgxn / bigX
    hgxp = hgxn[0] - 1
    on = [False] * n
    for i in range(n):
        if off[i]:
            on[i] = False
            continue
        on[i] = (hgxn[i] > hgxp)
        if on[i]:
            hgxp = hgxn[i]
    ittot = 1
    fit_gtspline(n, xs, yn, on, qt, jump, yat, en, FF)
    ena = en
    if FF:
        print('In slalom_tspline; failure flag raised in call to fit_gtspline')
        print('at initialization of A loop')
        return
    for ita in range(nita):
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
        for itb in range(nitb):
            fit_gtspline(n, xs, yn, on, qt, jump, yat, en, FF)
            if FF:
                print('In slalom_tspline; failure flag raised in call to fit_gtspline')
                print('at B loop, iterations ita, itb =', ita, itb)
                return
            ittot += 1
            fit_tspline()
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
            print('In slalom_tspline; exceeding the allocation of B iterations')
            return
        q = qt
        ya = yat
        if en >= ena:
            print('In slalom_tspline; energy failed to decrease')
            break
        ena = en
    if ita > nita:
        FF = True
        print('In slalom_tspline; exceeding the allocation of A iterations')
        return



