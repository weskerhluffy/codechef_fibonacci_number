//
//  main.cpp
//  shame
//
//  Created by ernesto alvarado on 05/09/19.
//  Copyright © 2019 ernesto alvarado. All rights reserved.
//

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <map>
using namespace std;

namespace QuadraticResidue
{
    typedef long long Long;
    
    Long powMod(Long a, Long n, Long m) {
        Long result = 1 % m;
        while (n > 0) {
            if ((n & 1) == 1) {
                result = result * a % m;
            }
            a = a * a % m;
            n >>= 1;
        }
        return result;
    }
    
    Long legendre(Long a, Long p) {
        return powMod(a, (p - 1) >> 1, p);
    }
    
    Long mod(Long a, Long m) {
        a %= m;
        if (a < 0) {
            a += m;
        }
        return a;
    }
    
    struct Combination {
        static Long p, omega;
        
        Long a, b;
        
        Combination(Long a, Long b): a(a % p), b(b % p) {}
    };
    
    Combination operator *(const Combination &p, const Combination &q) {
        int m = Combination::p;
        return Combination(p.a * q.a + p.b * q.b % m * Combination::omega,
                           p.a * q.b + q.a * p.b);
    }
    
    Combination powMod(Combination a, Long n) {
        Combination result(1, 0);
        while (n > 0) {
            if ((n & 1) == 1) {
                result = result * a;
            }
            a = a * a;
            n >>= 1;
        }
        return result;
    }
    
    Long Combination::p, Combination::omega;
    
    Long solve(Long a, Long p) {
        if (p == 2) {
            return 1;
        }
        if (legendre(a, p) + 1 == p) {
            return -1;
        }
        if ((((p + 1) >> 1) & 1) == 0) {
            return powMod(a, (p + 1) >> 2, p);
        }
        Long a_0 = -1;
        while (true) {
            a_0 = rand() % p;
            if (legendre(mod(a_0 * a_0 - a, p), p) + 1 == p) {
                break;
            }
        }
        Combination::p = p;
        Combination::omega = mod(a_0 * a_0 - a, p);
        //printf("%lld\n", Combination::omega);
        Combination ret = powMod(Combination(a_0, 1), (p + 1) >> 1);
        assert(ret.b == 0);
        return ret.a;
    }
    
    int getRoot(int a, int p) {
        a %= p;
        int x_0 = solve(a, p);
        if (x_0 == -1) {
            return -1;
        } else {
            return x_0;
        }
    }
};

inline void multiply(int &a, int b, int p)
{
    a = (long long)a * b % p;
}

inline int powmod(int x, int t, int mod)
{
    if (!t) {
        return 1 % mod;
    }
    int y = powmod(x, t >> 1, mod);
    y = (long long)y * y % mod;
    if (t & 1) {
        y = (long long)y * x % mod;
    }
    return y;
}

inline int reverse(int x, int p)
{
    return powmod(x, p - 2, p);
}

map<int, int> hash[2];

inline int solve(int a, int c, int p, int n, int sqrtP)
{
    int cDiv2 = (long long)c * reverse(2, p) % p;
    int temp = (long long)cDiv2 * cDiv2 % p;
    if (n & 1) {
        temp = (temp - 1 + p) % p;
    } else {
        temp = (temp + 1 + p) % p;
    }
    //fprintf(stderr, "	a = %d, c = %d, p = %d, cDiv2 = %d, temp = %d\n", a, c, p, cDiv2, temp);
    
    int root = QuadraticResidue::getRoot(temp, p);
    if (root == -1) {
        return -1;
    }
    //fprintf(stderr, "	root = %d\n", root);
    //a ^ n - c / 2 = root
    int bigA = powmod(a, sqrtP, p);
    int ret = -1;
    for (int iter = 0; iter < 2; ++ iter) {
        int target = (root + cDiv2) % p;
        //fprintf(stderr, "		target = %d\n", target);
        //bigA ^ x * a ^ y = target    mod p
        int mul = 1;
        for (int x = 0; x < sqrtP; ++ x) {
            //mul * a ^ y = target
            //a ^ y = target / mul
            int ay = (long long)target * reverse(mul, p) % p;
            int odd = n ^ (x * sqrtP & 1);
            //fprintf(stderr, "			odd = %d, ay = %d, target = %d, mul = %d\n", odd, ay, target, reverse(mul, p));
            if (hash[odd].count(ay)) {
                int candidate = x * sqrtP + hash[odd][ay];
                //fprintf(stderr, "			candidate = %d\n", candidate);
                if (ret == -1 || ret > candidate) {
                    ret = candidate;
                }
            }
            multiply(mul, bigA, p);
        }
        root = (p - root) % p;
    }
    return ret;
}

inline int solve(int c, int p)
{
    if (c == 0) {
        return 0;
    }
    if (c == 1) {
        return 1;
    }
    //fprintf(stderr, "solve %d %d\n", c, p);
    //x ^ 2 = 5 mod p
    int x = QuadraticResidue::getRoot(5, p); // sqrt(5)
    //fprintf(stderr, "sqrt(5) mod p = %d, assert %d\n", x, (long long)x * x % p);
    assert(x != -1);
    //a = (1 + sqrt(5)) / 2
    int a = (long long)(1 + x) * reverse(2, p) % p;
    //fib(n) = 1 / sqrt(5) (a ^ n - (-1 / a) ^ n)
    c = (long long)c * x % p;
    //fprintf(stderr, "newC = %d\n", c);
    int sqrtP = (int)sqrt(p) + 1;
    int mul = 1;
    hash[0].clear();
    hash[1].clear();
    for (int x = 0; x < sqrtP; ++ x) {
        if (!hash[x & 1].count(mul)) {
            hash[x & 1][mul] = x;
            //fprintf(stderr, "cached: %d %d: %d\n", x & 1, mul, x);
        }
        multiply(mul, a, p);
    }
    
    int ret = -1;
    for (int odd = 0; odd < 2; ++ odd) {
        int t = solve(a, c, p, odd, sqrtP);
        if (t != -1) {
            if (ret == -1 || ret > t) {
                ret = t;
            }
        }
    }
    return ret;
}

int main()
{
    freopen("input.txt", "r", stdin);
    
    int T;
    for (scanf("%d", &T); T --;) {
        int c, p;
        scanf("%d%d", &c, &p);
        printf("%d\n", solve(c, p));
    }
    return 0;
}
