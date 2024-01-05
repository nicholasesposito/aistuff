def eval_usplined(n, xs, p, q, x):
    if x <= xs[0]:
        xr = x - xs[0]
        y = p[0] + q[0] * xr
        dydx = q[0]
        return y, dydx
    if x >= xs[n-1]:
        xr = x - xs[n-1]
        y = p[n-1] + q[n-1] * xr
        dydx = q[n-1]
        return y, dydx
    for ib in range(1, n):
        if xs[ib] <= xs[ib-1]:
            continue
        if xs[ib] >= x:
            break
    ia = ib - 1
    xh = (xs[ib] - xs[ia]) / 2
    xr = x - xs[ia] - xh
    pm = (p[ia] + p[ib]) / 2
    qm = (p[ib] - p[ia]) / (xh * 2)
    qah = q[ia] / 2
    qbh = q[ib] / 2
    qxh = qah + qbh - qm
    qdh = qbh - qah
    shh = xh
    chh = 1
    sh = xr
    ch = 1
    shm = xr**3 / 6
    chm = xr**2 / 2
    shhm = xh**3 / 6
    chhm = xh**2 / 2
    xcmsh = xh**3 / 3
    qdh = qdh / shh
    qxh = qxh / xcmsh
    y = pm + xr * qm + qdh * (chm - chhm) + qxh * (xh * shm - xr * shhm)
    dydx = qm + qdh * sh + qxh * (xh * chm - shhm)
    return y, dydx

def eval_usplinedd(n, xs, p, q, x):
    if x <= xs[0]:
        xr = x - xs[0]
        y = p[0] + q[0] * xr
        dydx = q[0]
        ddydxx = 0
        return y, dydx, ddydxx
    if x >= xs[n-1]:
        xr = x - xs[n-1]
        y = p[n-1] + q[n-1] * xr
        dydx = q[n-1]
        ddydxx = 0
        return y, dydx, ddydxx
    for ib in range(1, n):
        if xs[ib] <= xs[ib-1]:
            continue
        if xs[ib] >= x:
            break
    ia = ib - 1
    xh = (xs[ib] - xs[ia]) / 2
    xr = x - xs[ia] - xh
    pm = (p[ia] + p[ib]) / 2
    qm = (p[ib] - p[ia]) / (xh * 2)
    qah = q[ia] / 2
    qbh = q[ib] / 2
    qxh = qah + qbh - qm
    qdh = qbh - qah
    shh = xh
    chh = 1
    sh = xr
    ch = 1
    shm = xr**3 / 6
    chm = xr**2 / 2
    shhm = xh**3 / 6
    chhm = xh**2 / 2
    xcmsh = xh**3 / 3
    qdh = qdh / shh
    qxh = qxh / xcmsh
    y = pm + xr * qm + qdh * (chm - chhm) + qxh * (xh * shm - xr * shhm)
    dydx = qm + qdh * sh + qxh * (xh * chm - shhm)
    ddydxx = qdh + qxh * xh * sh
    return y, dydx, ddydxx

def eval_usplineddd(n, xs, p, q, x):
    if x <= xs[0]:
        xr = x - xs[0]
        y = p[0] + q[0] * xr
        dydx = q[0]
        ddydxx = 0
        dddydxxx = 0
        return y, dydx, ddydxx, dddydxxx
    if x >= xs[n-1]:
        xr = x - xs[n-1]
        y = p[n-1] + q[n-1] * xr
        dydx = q[n-1]
        ddydxx = 0
        dddydxxx = 0
        return y, dydx, ddydxx, dddydxxx
    for ib in range(1, n):
        if xs[ib] <= xs[ib-1]:
            continue
        if xs[ib] >= x:
            break
    ia = ib - 1
    xh = (xs[ib] - xs[ia]) / 2
    xr = x - xs[ia] - xh
    pm = (p[ia] + p[ib]) / 2
    qm = (p[ib] - p[ia]) / (xh * 2)
    qah = q[ia] / 2
    qbh = q[ib] / 2
    qxh = qah + qbh - qm
    qdh = qbh - qah
    shh = xh
    chh = 1
    sh = xr
    ch = 1
    shm = xr**3 / 6
    chm = xr**2 / 2
    shhm = xh**3 / 6
    chhm = xh**2 / 2
    xcmsh = xh**3 / 3
    qdh = qdh / shh
    qxh = qxh / xcmsh
    y = pm + xr * qm + qdh * (chm - chhm) + qxh * (xh * shm - xr * shhm)
    dydx = qm + qdh * sh + qxh * (xh * chm - shhm)
    ddydxx = qdh + qxh * xh * sh
    dddydxxx = qxh * xh
    return y, dydx, ddydxx, dddydxxx

def eval_iuspline(n, xs, p, q, m, x):
    if x <= xs[0]:
        xr = x - xs[0]
        y = p[0] * xr + q[0] * xr**2 / 2
        return y
    if x >= xs[n-1]:
        xr = x - xs[n-1]
        y = m[n-1] + p[n-1] * xr + q[n-1] * xr**2 / 2
        return y
    for ib in range(1, n):
        if xs[ib] <= xs[ib-1]:
            continue
        if xs[ib] >= x:
            break
    ia = ib - 1
    xh = (xs[ib] - xs[ia]) / 2
    xr = x - xs[ia] - xh
    t2 = xh**2 / 2
    t3 = t2 * xh / 3
    pa = (p[ib] + p[ia]) / 2
    pd = (p[ib] - p[ia]) / xh
    qa = (q[ib] + q[ia]) / 2
    qd = (q[ib] - q[ia]) / xh
    c = qd
    a = pa - c * t2
    d = (qa - pd) * 3/2 / t2
    b = qa - d * t2
    t2 = xr**2 / 2
    t3 = t2 * xr / 3
    t4 = t3 * xr / 4
    y = m[ia] + a * xr + b * t2 + c * t3 + d * t4
    return y

def best_tslalom(nh, mh, doru, hgts, hs, halfgate, bigT):
    hgtp = 0
    hgtn = [[0] * mh for _ in range(2)]
    hn = [[[0] * mh for _ in range(2)] for _ in range(2)]
    code = [0] * mh
    mode = [0] * mh
    bend = [0] * (mh * 2)
    off = [False] * (mh * 2)
    FF = False
    set_gates(nh, mh, doru, hgts, hs, hgtn, hn, code, FF)
    if FF:
        print('In best_tslalom; failure flag was raised in call to set_gates')
        return
    count_routes(mh, code, route_count, FF)
    maxrts = max(maxrts, route_count)
    if FF:
        print('In best_tslalom; failure flag was raised in call to count_routes')
        return
    if route_count > 4:
        list_routes(mh, code)
    enbest = float('inf')
    flag = True
    for k in range(1, i+1):
        next_route(mh, code, mode, flag)
        if flag:
            flag = False
            break
        set_posts(mh, mode, hgtn, hn, bend, hgtp, hp, off)
        slalom_tspline(m, bend, hgtp, hp, off, bigT, q, ya, en, ita, maxitb, ittot, FF)
        en = en / enbase_t(tspan, hspan)
        maxita = max(maxita, ita)
        maxit = max(maxit, ittot)
        if FF:
            print('In best_tslalom; failure flag was raised in call to slalom_tspline')
            return
        if en < enbest:
            modebest = mode
            enbest = en
            qbest = q
            yabest = ya

def best_uslalom(nh, mh, doru, hgts, hs, halfgate):
    hgtp = 0
    hgtn = [[0] * mh for _ in range(2)]
    hn = [[[0] * mh for _ in range(2)] for _ in range(2)]
    code = [0] * mh
    mode = [0] * mh
    bend = [0] * (mh * 2)
    off = [False] * (mh * 2)
    FF = False
    set_gates(nh, mh, doru, hgts, hs, hgtn, hn, code, FF)
    if FF:
        print('In best_uslalom; failure flag was raised in call to set_gates')
        return
    count_routes(mh, code, route_count, FF)
    maxrts = max(maxrts, route_count)
    if FF:
        print('In best_uslalom; failure flag was raised in call to count_routes')
        return
    if route_count > 4:
        list_routes(mh, code)
    enbest = float('inf')
    flag = True
    for k in range(1, i+1):
        next_route(mh, code, mode, flag)
        if flag:
            flag = False
            break
        set_posts(mh, mode, hgtn, hn, bend, hgtp, hp, off)
        slalom_uspline(m, bend, hgtp, hp, off, halfgate, q, ya, en, ita, maxitb, ittot, FF)
        maxita = max(maxita, ita)
        maxit = max(maxit, ittot)
        if FF:
            print('In best_uslalom; failure flag was raised in call to slalom_uspline')
            return
        if en < enbest:
            modebest = mode
            enbest = en
            qbest = q
            yabest = ya

def count_gates(nh, hgts):
    hgtp = hgts[0] - 1
    mh = 0
    for i in range(nh):
        if hgts[i] <= hgtp:
            continue
        mh += 1
        hgtp = hgts[i]
    return mh



