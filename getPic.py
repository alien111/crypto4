import random
import time
import matplotlib.pyplot as plt
import numpy as np

A = random.randint(1000000000, 10000000000)
B = random.randint(1000000000, 10000000000)

def ellipticCurve(x, y, p):
    return (y ** 2) % p == (x ** 3 + (A % p) * x + (B % p)) % p

def searchPoints():
    points = []
    for x in range(p):
        for y in range(p):
            if ellipticCurve(x, y, p):
                points.append((x, y))
    return points

def extendedEuclideanAlgorithm(a, b):
    s, old_s = 0, 1
    t, old_t = 1, 0
    r, old_r = b, a

    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t

    return old_r, old_s, old_t


def inverse(n, p):
    gcd, x, y = extendedEuclideanAlgorithm(n, p)
    assert (n * x + p * y) % p == gcd

    if gcd != 1:
        raise ValueError(
            '{} has no multiplicative inverse '
            'modulo {}'.format(n, p))
    else:
        return x % p

def addPoints(p1, p2, p):
    x1, y1 = p1[0], p1[1]
    x2, y2 = p2[0], p2[1]

    if p1 == (0, 0):
        return p2
    elif p2 == (0, 0):
        return p1
    elif x1 == x2 and y1 != y2:
        return (0, 0)

    

    if p1 == p2:
        m = ((3 * x1 ** 2 + (A % p)) * inverse(2 * y1, p)) % p
    else:
        m = ((y1 - y2) * inverse(x1 - x2, p)) % p
        

    x3 = (m ** 2 - x1 - x2) % p
    y3 = (y1 + m * (x3 - x1)) % p

    return [x3, -y3 % p]


def pointOrder(point, p):
    i = 1
    new_point = addPoints(point, point, p)
    while new_point != (0, 0):
        new_point = addPoints(new_point, point, p)
        i += 1

    return i

def sieve(n):
    primes = 2 * [False] + (n - 1) * [True]
    for i in range(2, int(n ** 0.5 + 1.5)):
        for j in range(i * i, n + 1, i):
            primes[j] = False
    return [prime for prime, checked in enumerate(primes) if checked]

def throughPrimes(p):
    start = time.time()
    points = searchPoints()
    points_num = len(points)
    point = random.choice(points)
    order = pointOrder(point, p)
    stopTime = time.time() - start
    return stopTime

primes = sieve(4000)
stat = []

for i in range(2, len(primes)):
    p = primes[i]
    #print(p)
    stat.append(throughPrimes(p))

plt.scatter(primes[2:], stat)
plt.show()
