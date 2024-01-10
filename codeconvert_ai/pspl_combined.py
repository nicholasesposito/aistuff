import numpy as np

halfgate = np.float(30.0)
bigT = np.float(120.0)
heps = np.float(0.01)

def expm(x):
    x = np.float(x)
    if abs(x) > 0.5:
        e = np.exp(x) - 1.0
    else:
        p = np.float(x)
        e = np.float(p)
        for i in range(2, 20):
            p = (p * x) / i
            e = e + p
            if abs(p) <= abs(e * np.finfo(float).eps):
                return np.float(e)

def expmm(x):
    x = np.float(x)
    if abs(x) > 0.5:
        e = np.exp(x) - 1.0 - x
    else:
        p = np.float(x * x * 0.5)
        e = np.float(p)
        for i in range(3, 26):
            p = (p * x) / i
            e = e + p
            if abs(p) <= abs(e * np.finfo(float).eps):
                return np.float(e)

def coshm(x):
    x = np.float(x)
    return np.float(2) * np.sinh(x * 0.5) ** 2

def sinhm(x):
    x = np.float(x)
    if abs(x) > 0.5:
        s = np.sinh(x) - x
    else:
        p = np.float(x ** 3 / 6)
        s = np.float(p)
        xx = np.float(x * x)
        for i in range(5, 20, 2):
            p = p * xx / (i * (i - 1.0))
            s = s + p
            if abs(p) <= abs(s * np.finfo(float).eps):
                return np.float(s)

def coshmm(x):
    x = np.float(x)
    xh = x * 0.5
    return sinhm(xh) * (2 * np.sinh(xh) + x)

def xcms(x):
    x = np.float(x)
    if abs(x) > 0.5:
        e = np.float(x * coshm(x) - sinhm(x))
    else:
        p = np.float(x ** 3 / 3)
        e = np.float(p)
        xx = x * x
        for i in range(2, 16):
            i2 = i * 2
            p = p * xx / (i2 * (i2 + 1))
            e = e + i * p
            if abs(p) <= abs(e * np.finfo(float).eps):
                return np.float(e)

def enbase_t(tspan, hspan):
    tspan = np.float(tspan)
    hspan = np.float(hspan)
    if tspan < 0.0:
        raise ValueError("In enbase_t; thspan must be positive")
    if hspan == 0.0:
        return 1.0
    return (hspan ** 2 / expmm(-tspan)) * 0.5


def tbnewton(nh, m, hgts, hs, hgtp, p, q):
    nit = int(12)

    nh = int(nh)
    m = int(m)
    hgts = np.float(hgts)
    hs = np.float(hs)
    hgtp = np.float(hgtp)
    p = np.float(p) 
    q = np.float(q) 
    gate = 2 * halfgate / bigT
    tr = hgtp * halfgate / bigT
    for i in range(nh):
        tee = hgts[i] * halfgate / bigT
        he = hs[i]
        it = 1
        while it <= 12:
            eval_tspline(m, tr, p, q, tee, hac, dhadt)
            if it == 1:
                dhdt[i] = dhadt / bigT
            if dhadt == 0:
                break
            dh = hac - he
            dt = -dh / dhadt
            if abs(dt) > gate:
                print("WARNING! In tbnewton; i,it,dt/gate =", i, it, dt / gate)
                break
            if abs(dh) < np.finfo(float).eps:
                dhdt[i] = dhadt / bigT
                break
            tee = tee + dt
            it = it + 1
        FF = it > 12
        if FF:
            print("In tbnewton; Newton iterations seem not to be converging at i=", i)
            print("tee,he,hac,heps,dhadt:", tee, he, hac, np.finfo(float).eps, dhadt)
        te[i] = tee

        return te, dhdt, FF

def ubnewton(nh, m, hgts, hs, hgtp, p, q):
    nit = int(12)

    nh = int(nh)
    m = int(m)
    hgts = np.float(hgts)
    hs = np.float(hs)
    hgtp = np.float(hgtp)
    p = np.float(p)
    q = np.float(q)

    gate = 2 * halfgate
    tr = hgtp * halfgate
    for i in range(nh):
        tee = hgts[i] * halfgate
        he = hs[i]
        it = 1
        while it <= 12:
            eval_uspline(m, tr, p, q, tee, hac, dhadt)  #NE CHECK
            if it == 1:
                dhdt[i] = dhadt
            if dhadt == 0:
                break
            dh = hac - he
            dt = -dh / dhadt
            if abs(dt) > gate:
                print("WARNING! In ubnewton; i,it,dt/gate =", i, it, dt / gate)
                break
            if abs(dh) < np.finfo(float).eps:
                dhdt[i] = dhadt
                break
            tee = tee + dt
            it = it + 1
        FF = it > 12
        if FF:
            print("In ubnewton; Newton iterations seem not to be converging at i=", i)
            print("tee,he,hac,heps,dhadt:", tee, he, hac, np.finfo(float).eps, dhadt)
        te[i] = tee

    return te, dhdt, FF

def fit_gtspline(n, xs, ys, on):
    m = 0
    xa = np.zeros(n)
    ya = np.zeros(n)
    qa = np.zeros(n)
    ja = np.zeros(n)
    for i in range(n):
        if on[i]:
            m = m + 1
            xa[m-1] = xs[i]
            ya[m-1] = ys[i]
    fit_tspline(m, xa[:m], ya[:m], qa[:m], ja[:m], en, FF)
    if FF:
        print("In fit_gtspline; failure flag raised at call to fit_tspline")
        return
    k = 0
    for i in range(n):
        if on[i]:
            k = k + 1
            q[i] = qa[k-1]
            j[i] = ja[k-1]
            yac[i] = ys[i]
        else:
            eval_tsplined(m, xa[:m], ya[:m], qa[:m], xs[i], yac[i], q[i])
            j[i] = 0

    return q, j, yac, en, FF

def fit_tspline(n, xs, p, q, j, en, FF):
    if n < 1:
        raise ValueError("In fit_tspline; size of data array must be positive")
    if n == 1:
        q[:] = 0
        j[:] = 0
        en = 0
        return
    for i in range(1, n):
        if xs[i-1] >= xs[i]:
            FF = True
            print("In fit_tspline; xs data must increase strictly monotonically")
            return
    qq = np.zeros((n, 2))
    difp = np.zeros(n-1)
    cpp = np.zeros(n-1)
    cqp = np.zeros(n-1)
    sumq = np.zeros(n-1)
    for i in range(n-1):
        ip = i + 1
        difp[i] = p[ip] - p[i]
        x = (xs[ip] - xs[i]) * 0.5
        ch = np.cosh(x)
        sh = np.sinh(x)
        xcmsx2 = xcms(x) * 2
        egg = x * sh / xcmsx2
        ehh = ch / (2 * sh)
        ccc = egg + ehh
        cpp[i] = ch / xcmsx2
        cqp[i] = -difp[i] * sh / xcmsx2
        qq[i, 0] = qq[i, 0] + ccc
        qq[ip, -1] = qq[ip, -1] + egg - ehh
        qq[ip, 0] = qq[ip, 0] + ccc
    qq[0, 0] = qq[0, 0] + 1
    qq[n-1, 0] = qq[n-1, 0] + 1
    q[:n-1] = -cqp
    q[n-1] = 0
    q[1:n] = q[1:n] - cqp
    ldltb(n, 1, qq)
    ltdlbv(n, 1, qq, q)
    sumq[:] = q[:n-1] + q[1:n]
    en = 0.5 * (np.dot(difp**2, cpp) + np.dot(sumq, cqp))
    sb = q[0]
    for i in range(n-1):
        ip = i + 1
        x = (xs[ip] - xs[i]) * 0.5
        xcmsx2 = xcms(x) * 2
        ch = np.cosh(x)
        sh = np.sinh(x)
        sap = (sh * sumq[i] - ch * difp[i]) / xcmsx2
        sa = sap + q[i]
        j[i] = sa - sb
        sb = sap + q[ip]
    j[n-1] = q[n-1] - sb

def int_tspline(n, xs, p, q, m):
    m[0] = 0
    e = 0
    for i in range(n-1):
        ip = i + 1
        x = (xs[ip] - xs[i]) * 0.5
        t2 = x * x * 0.5
        shx = np.sinh(x)
        chmx = coshm(x)
        shmx = sinhm(x)
        chmmx = coshmm(x)
        xcmsx = xcms(x)
        pa = (p[ip] + p[i]) * 0.5
        pd = (p[ip] - p[i]) * 0.5 / x
        qa = (q[ip] + q[i]) * 0.5
        qd = (q[ip] - q[i]) * 0.5 / shx
        c = qd
        a = pa - c * chmx
        d = (qa - pd) * x / xcmsx
        b = qa - d * chmx
        m[i] = e + a * x - b * t2 + c * shmx - d * chmmx
        e = e + 2 * (a * x + c * shmx)
    m[n-1] = e

def fit_guspline(n, xs, ys, on, q, j, yac, en, FF):
    m = 0
    xa = np.zeros(n)
    ya = np.zeros(n)
    qa = np.zeros(n)
    ja = np.zeros(n)
    for i in range(n):
        if on[i]:
            m = m + 1
            xa[m-1] = xs[i]
            ya[m-1] = ys[i]
    fit_uspline(m, xa[:m], ya[:m], qa[:m], ja[:m], en, FF)
    if FF:
        print("In fit_guspline; failure flag raised at call to fit_uspline")
        return
    k = 0
    for i in range(n):
        if on[i]:
            k = k + 1
            q[i] = qa[k-1]
            j[i] = ja[k-1]
            yac[i] = ys[i]
        else:
            eval_usplined(m, xa[:m], ya[:m], qa[:m], xs[i], yac[i], q[i])
            j[i] = 0

def fit_uspline(n, xs, p, q, j, en, FF):
    if n < 1:
        raise ValueError("In fit_uspline; size of data array must be positive")
    if n == 1:
        q[:] = 0
        j[:] = 0
        en = 0
        return
    for i in range(1, n):
        if xs[i-1] >= xs[i]:
            FF = True
            print("In fit_uspline; xs data must increase strictly monotonically")
            return
    qq = np.zeros((n, 2))
    difp = np.zeros(n-1)
    cpp = np.zeros(n-1)
    cqp = np.zeros(n-1)
    sumq = np.zeros(n-1)
    for i in range(n-1):
        ip = i + 1
        difp[i] = p[ip] - p[i]
        x = (xs[ip] - xs[i]) * 0.5
        ch = np.cosh(x)
        sh = np.sinh(x)
        xcmsx2 = xcms(x) * 2
        egg = x * sh / xcmsx2
        ehh = ch / (2 * sh)
        ccc = egg + ehh
        cpp[i] = ch / xcmsx2
        cqp[i] = -difp[i] * sh / xcmsx2
        qq[i, 0] = qq[i, 0] + ccc
        qq[ip, -1] = qq[ip, -1] + egg - ehh
        qq[ip, 0] = qq[ip, 0] + ccc
    qq[0, 0] = qq[0, 0] + 1
    qq[n-1, 0] = qq[n-1, 0] + 1
    q[:n-1] = -cqp
    q[n-1] = 0
    q[1:n] = q[1:n] - cqp
    ldltb(n, 1, qq)
    ltdlbv(n, 1, qq, q)
    sumq[:] = q[:n-1] + q[1:n]
    en = 0.5 * (np.dot(difp**2, cpp) + np.dot(sumq, cqp))
    sb = q[0]
    for i in range(n-1):
        ip = i + 1
        x = (xs[ip] - xs[i]) * 0.5
        xcmsx2 = xcms(x) * 2
        ch = np.cosh(x)
        sh = np.sinh(x)
        sap = (sh * sumq[i] - ch * difp[i]) / xcmsx2
        sa = sap + q[i]
        j[i] = sa - sb
        sb = sap + q[ip]
    j[n-1] = q[n-1] - sb

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
