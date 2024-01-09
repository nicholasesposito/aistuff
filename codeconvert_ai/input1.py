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



