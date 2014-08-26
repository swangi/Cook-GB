import numpy as np
from fractions import gcd
shif = 1
anga_100 = []
anga_110 = []
angle_100 = []
angle_110 = []
mlist_100 = []
nlist_100 = []
sigmalist_100 = []
mlist_110 = []
nlist_110 = []
sigmalist_110 = []
tol_100 = 0.735
tol_110 = 1.4
# 100 tilt #####################
for m in np.arange(1, 30):
    for n in np.arange(1, 30):
        if n < m:
            continue
        pnp = 0
        if gcd(m, n) != 1:
            continue
        ang = 2 * np.degrees(np.arctan2(m, n))
        anga_100.append(ang)
        for a in anga_100:
            if (ang - a) != 0 and abs(ang - a) < tol_100:
                pnp = 1
        if pnp == 1:
            continue
        Sigma = m ** 2 + n ** 2 + 0 ** 2
        while Sigma % 2 == 0:
            Sigma = Sigma / 2
        mlist_100.append(m)
        nlist_100.append(n)
        angle = 2 * np.degrees(np.arctan2(m, n))
        angle_100.append(ang)
        sigmalist_100.append(Sigma)
# 110 tilt #####################
for m in np.arange(1, 20):
    for n in np.arange(1, 20):
        pnp = 0
        if gcd(m, n) != 1:
            continue
        ang = 2 * np.degrees(np.arctan2(np.sqrt(2) * m, n))
        anga_110.append(ang)
        for a in anga_110:
            if (ang - a) != 0 and abs(ang - a) < tol_110:
                pnp = 1
        if pnp == 1:
            continue
        Sigma = m ** 2 + m ** 2 + n ** 2
        while Sigma % 2 == 0:
            Sigma = Sigma / 2
        mlist_110.append(m)
        nlist_110.append(n)
        angle = 2 * np.degrees(np.arctan2(m, n))
        angle_110.append(ang)
        sigmalist_110.append(Sigma)
if __name__ == '__main__':
    for i, m in enumerate(mlist_110):
        print "Sigma%d[%d, %d, %d]\t%.5f\n" % (sigmalist_100[i], mlist_100[i], mlist_100[i], nlist_100[i], angle_100[i])
    # print mlist_110
    # print nlist_110
    # print anga_110
