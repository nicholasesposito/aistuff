def fit_uspline(n, xs, p):
    import numpy as np
    q = np.zeros(n)
    j = np.zeros(n)
    en = 0
    FF = False
    if n < 1:
        raise ValueError('In fit_uspline; size of data array must be positive')
    if n == 1:
        q = np.zeros(n)
        j = np.zeros(n)
        en = 0
        return q, j, en, FF
    for i in range(1, n):
        if xs[i-1] >= xs[i]:
            FF = True
            print('In fit_uspline; xs data must increase strictly monotonically')
            return q, j, en, FF
    qq = np.zeros((n, 2))
    for i in range(n-1):
        ip = i + 1
        difp = p[ip] - p[i]
        x2 = xs[ip] - xs[i]
        x = x2 / 2
        xcmsx2 = (x ** 3) * 2 / 3
        ccc = 2 / x
        cpp = 1 / xcmsx2
        cqp = -difp * x / xcmsx2
        qq[i, 0] += ccc
        qq[ip, -1] += 1 / x
        qq[ip, 0] += ccc
        q[0:n-1] = -cqp
        q[n-1] = 0
        q[1:n] -= cqp
    ldltb(n, 1, qq)
    ltdlbv(n, 1, qq, q)
    sumq = q[0:n-1] + q[1:n]
    en = (np.dot(difp ** 2, cpp) + np.dot(sumq, cqp)) / 2
    sb = 0
    for i in range(n-1):
        ip = i + 1
        x = (xs[ip] - xs[i]) / 2
        xcmsx2 = (x ** 3) * 2 / 3
        sa = (x * sumq[i] - difp[i]) / xcmsx2
        j[i] = sa - sb
        sb = sa
    j[n-1] = -sb
    return q, j, en, FF

def int_uspline(n, xs, p, q):
    import numpy as np
    m = np.zeros(n)
    u3o2 = 3 / 2
    e = 0
    for i in range(n-1):
        ip = i + 1
        x = (xs[ip] - xs[i]) / 2
        t2 = x ** 2 / 2
        t3 = t2 * x / 3
        t4 = t3 * x / 4
        pa = (p[ip] + p[i]) / 2
        pd = (p[ip] - p[i]) / (2 * x)
        qa = (q[ip] + q[i]) / 2
        qd = (q[ip] - q[i]) / (2 * x)
        c = qd
        a = pa - c * t2
        d = (qa - pd) * u3o2 / t2
        b = qa - d * t2
        m[i] = e + a * x - b * t2 + c * t3 - d * t4
        e = e + 2 * (a * x + c * t3)
    m[n-1] = e
    return m

def eval_tspline(n, xs, p, q, x):
    import numpy as np
    ia = 0
    ib = 0
    xr = 0
    xh = 0
    pm = 0
    qm = 0
    qah = 0
    qbh = 0
    qxh = 0
    qdh = 0
    shh = 0
    chh = 0
    sh = 0
    ch = 0
    xcmsh = 0
    shm = 0
    chm = 0
    shhm = 0
    chhm = 0
    if x <= xs[0]:
        xr = x - xs[0]
        y = p[0] + q[0] * np.exp(xr)
        return y
    if x >= xs[n-1]:
        xr = x - xs[n-1]
        y = p[n-1] - q[n-1] * np.exp(-xr)
        return y
    for ib in range(1, n):
        if xs[ib] <= xs[ib-1]:
            continue
        if xs[ib] >= x:
            break
    ia = ib - 1
    xh = (xs[ib] - xs[ia]) / 2
    xr = x - xs[ia] - xh
    pm = (p[ib] + p[ia]) / 2
    qm = (p[ib] - p[ia]) / (2 * xh)
    qah = q[ia] / 2
    qbh = q[ib] / 2
    qxh = qah + qbh - qm
    qdh = qbh - qah
    shh = np.sinh(xh)
    chh = np.cosh(xh)
    sh = np.sinh(xr)
    ch = np.cosh(xr)
    shm = np.sinh(xr) / 6
    chm = np.cosh(xr) / 2
    shhm = np.sinh(xh) / 6
    chhm = np.cosh(xh) / 2
    xcmsh = (xh ** 3) / 3
    qdh = qdh / shh
    qxh = qxh / xcmsh
    y = pm + xr * qm + qdh * (chm - chhm) + qxh * (xh * shm - xr * shhm)
    return y

def eval_tsplined(n, xs, p, q, x):
    import numpy as np
    ia = 0
    ib = 0
    xr = 0
    xh = 0
    pm = 0
    qm = 0
    qah = 0
    qbh = 0
    qxh = 0
    qdh = 0
    shh = 0
    chh = 0
    sh = 0
    ch = 0
    xcmsh = 0
    shm = 0
    chm = 0
    shhm = 0
    chhm = 0
    if x <= xs[0]:
        xr = x - xs[0]
        qemxr = q[0] * np.exp(xr)
        y = p[0] + qemxr
        dydx = qemxr + q[0]
        return y, dydx
    if x >= xs[n-1]:
        xr = x - xs[n-1]
        qemxr = q[n-1] * np.exp(-xr)
        y = p[n-1] - qemxr
        dydx = qemxr + q[n-1]
        return y, dydx
    for ib in range(1, n):
        if xs[ib] <= xs[ib-1]:
            continue
        if xs[ib] >= x:
            break
    ia = ib - 1
    xh = (xs[ib] - xs[ia]) / 2
    xr = x - xs[ia] - xh
    pm = (p[ib] + p[ia]) / 2
    qm = (p[ib] - p[ia]) / (2 * xh)
    qah = q[ia] / 2
    qbh = q[ib] / 2
    qxh = qah + qbh - qm
    qdh = qbh - qah
    shh = np.sinh(xh)
    chh = np.cosh(xh)
    sh = np.sinh(xr)
    ch = np.cosh(xr)
    shm = np.sinh(xr) / 6
    chm = np.cosh(xr) / 2
    shhm = np.sinh(xh) / 6
    chhm = np.cosh(xh) / 2
    xcmsh = (xh ** 3) / 3
    qdh = qdh / shh
    qxh = qxh / xcmsh
    y = pm + xr * qm + qdh * (chm - chhm) + qxh * (xh * shm - xr * shhm)
    dydx = qm + qdh * sh + qxh * (xh * chm - shhm)
    return y, dydx

def eval_tsplinedd(n, xs, p, q, x):
    import numpy as np
    ia = 0
    ib = 0
    xr = 0
    xh = 0
    pm = 0
    qm = 0
    qah = 0
    qbh = 0
    qxh = 0
    qdh = 0
    shh = 0
    chh = 0
    sh = 0
    ch = 0
    xcmsh = 0
    shm = 0
    chm = 0
    shhm = 0
    chhm = 0
    qemxr = 0
    if x <= xs[0]:
        xr = x - xs[0]
        qemxr = q[0] * np.exp(xr)
        y = p[0] + qemxr
        dydx = qemxr + q[0]
        ddydxx = dydx
        return y, dydx, ddydxx
    if x >= xs[n-1]:
        xr = x - xs[n-1]
        qemxr = q[n-1] * np.exp(-xr)
        y = p[n-1] - qemxr
        dydx = qemxr + q[n-1]
        ddydxx = -dydx
        return y, dydx, ddydxx
    for ib in range(1, n):
        if xs[ib] <= xs[ib-1]:
            continue
        if xs[ib] >= x:
            break
    ia = ib - 1
    xh = (xs[ib] - xs[ia]) / 2
    xr = x - xs[ia] - xh
    pm = (p[ib] + p[ia]) / 2
    qm = (p[ib] - p[ia]) / (2 * xh)
    qah = q[ia] / 2
    qbh = q[ib] / 2
    qxh = qah + qbh - qm
    qdh = qbh - qah
    shh = np.sinh(xh)
    chh = np.cosh(xh)
    sh = np.sinh(xr)
    ch = np.cosh(xr)
    shm = np.sinh(xr) / 6
    chm = np.cosh(xr) / 2
    shhm = np.sinh(xh) / 6
    chhm = np.cosh(xh) / 2
    xcmsh = (xh ** 3) / 3
    qdh = qdh / shh
    qxh = qxh / xcmsh
    y = pm + xr * qm + qdh * (chm - chhm) + qxh * (xh * shm - xr * shhm)
    dydx = qm + qdh * sh + qxh * (xh * chm - shhm)
    ddydxx = qdh * ch + qxh * xh * sh
    return y, dydx, ddydxx

def eval_tsplineddd(n, xs, p, q, x):
    import numpy as np
    ia = 0
    ib = 0
    xr = 0
    xh = 0
    pm = 0
    qm = 0
    qah = 0
    qbh = 0
    qxh = 0
    qdh = 0
    shh = 0
    chh = 0
    sh = 0
    ch = 0
    xcmsh = 0
    shm = 0
    chm = 0
    shhm = 0
    chhm = 0
    qemxr = 0
    if x <= xs[0]:
        xr = x - xs[0]
        qemxr = q[0] * np.exp(xr)
        y = p[0] + qemxr
        dydx = qemxr + q[0]
        ddydxx = dydx
        dddydxxx = dydx
        return y, dydx, ddydxx, dddydxxx
    if x >= xs[n-1]:
        xr = x - xs[n-1]
        qemxr = q[n-1] * np.exp(-xr)
        y = p[n-1] - qemxr
        dydx = qemxr + q[n-1]
        ddydxx = -dydx
        dddydxxx = dydx
        return y, dydx, ddydxx, dddydxxx
    for ib in range(1, n):
        if xs[ib] <= xs[ib-1]:
            continue
        if xs[ib] >= x:
            break
    ia = ib - 1
    xh = (xs[ib] - xs[ia]) / 2
    xr = x - xs[ia] - xh
    pm = (p[ib] + p[ia]) / 2
    qm = (p[ib] - p[ia]) / (2 * xh)
    qah = q[ia] / 2
    qbh = q[ib] / 2
    qxh = qah + qbh - qm
    qdh = qbh - qah
    shh = np.sinh(xh)
    chh = np.cosh(xh)
    sh = np.sinh(xr)
    ch = np.cosh(xr)
    shm = np.sinh(xr) / 6
    chm = np.cosh(xr) / 2
    shhm = np.sinh(xh) / 6
    chhm = np.cosh(xh) / 2
    xcmsh = (xh ** 3) / 3
    qdh = qdh / shh
    qxh = qxh / xcmsh
    y = pm + xr * qm + qdh * (chm - chhm) + qxh * (xh * shm - xr * shhm)
    dydx = qm + qdh * sh + qxh * (xh * chm - shhm)
    ddydxx = qdh * ch + qxh * xh * sh
    dddydxxx = qdh * sh + qxh * xh * ch
    return y, dydx, ddydxx, dddydxxx

def eval_itspline(n, xs, p, q, m, x):
    import numpy as np
    ia = 0
    ib = 0
    xr = 0
    xh = 0
    pm = 0
    qm = 0
    qah = 0
    qbh = 0
    qxh = 0
    qdh = 0
    shh = 0
    chh = 0
    sh = 0
    ch = 0
    xcmsh = 0
    shm = 0
    chm = 0
    shhm = 0
    chhm = 0
    if x <= xs[0]:
        xr = x - xs[0]
        y = p[0] * xr + q[0] * np.expmm(xr)
        return y
    if x >= xs[n-1]:
        xr = x - xs[n-1]
        y = m[n-1] + p[n-1] * xr + q[n-1] * np.expmm(-xr)
        return y
    for ib in range(1, n):
        if xs[ib] <= xs[ib-1]:
            continue
        if xs[ib] >= x:
            break
    ia = ib - 1
    xh = (xs[ib] - xs[ia]) / 2
    xr = x - xs[ia] - xh
    pm = (p[ib] + p[ia]) / 2
    qm = (p[ib] - p[ia]) / (2 * xh)
    qah = q[ia] / 2
    qbh = q[ib] / 2
    qxh = qah + qbh - qm
    qdh = qbh - qah
    shh = xh
    chh = 1
    sh = xr
    ch = 1
    shm = xr ** 3 / 6
    chm = xr ** 2 / 2
    shhm = xh ** 3 / 6
    chhm = xh ** 2 / 2
    xcmsh = xh ** 3 / 3
    qdh = qdh / shh
    qxh = qxh / xcmsh
    y = m[ia] + pm * xr + qxh * (xh * shm - xr * shhm)
    return y

def eval_uspline(n, xs, p, q, x):
    import numpy as np
    ia = 0
    ib = 0
    xr = 0
    xh = 0
    pm = 0
    qm = 0
    qah = 0
    qbh = 0
    qxh = 0
    qdh = 0
    shh = 0
    chh = 0
    sh = 0
    ch = 0
    xcmsh = 0
    shm = 0
    chm = 0
    shhm = 0
    chhm = 0
    if x <= xs[0]:
        xr = x - xs[0]
        y = p[0] + q[0] * xr
        return y
    if x >= xs[n-1]:
        xr = x - xs[n-1]
        y = p[n-1] + q[n-1] * xr
        return y
    for ib in range(1, n):
        if xs[ib] <= xs[ib-1]:
            continue
        if xs[ib] >= x:
            break
    ia = ib - 1
    xh = (xs[ib] - xs[ia]) / 2
    xr = x - xs[ia] - xh
    pm = (p[ib] + p[ia]) / 2
    qm = (p[ib] - p[ia]) / (2 * xh)
    qah = q[ia] / 2
    qbh = q[ib] / 2
    qxh = qah + qbh - qm
    qdh = qbh - qah
    shh = xh
    chh = 1
    sh = xr
    ch = 1
    xcmsh = xh ** 3 / 3
    shm = xr ** 3 / 6
    chm = xr ** 2 / 2
    shhm = xh ** 3 / 6
    chhm = xh ** 2 / 2
    qdh = qdh / shh
    qxh = qxh / xcmsh
    y = pm + xr * qm + qdh * (chm - chhm) + qxh * (xh * shm - xr * shhm)
    return y



