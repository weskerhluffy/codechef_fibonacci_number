//
//  main.cpp
//  shame1
//
//  Created by ernesto alvarado on 02/10/19.
//  Copyright © 2019 ernesto alvarado. All rights reserved.
//
#include <ctime>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdio.h>
#include <iostream>
#include <map>

#include <time.h>
#include <algorithm>
#define Rg register
#define ll long long
using namespace std;
const int inf=0x7fffffff;
ll w,a,mod;
inline int dec(int x,int y){return x<y?x-y+mod:x-y;}
inline int inc(int x,int y){return 0ll+x+y>=mod?0ll+x+y-mod:x+y;}
inline int mul(int x,int y){return 1ll*x*y-1ll*x*y/mod*mod;}
inline int qpow(Rg int x,Rg int p,Rg int s=1){
    for(;p;p>>=1,x=mul(x,x)) if(p&1) s=mul(s,x); return s;
}
struct cp{ int x,y; inline cp(Rg int xx,Rg int yy){x=xx,y=yy;}
    inline cp operator *(const cp& b)const{ //有点鬼畜不知原理的向量乘操作呢
        return cp(inc(mul(x,b.x),mul(w,mul(y,b.y))),inc(mul(x,b.y),mul(y,b.x)));
    }
};
inline int qpow(Rg cp x,Rg int p){ Rg cp s(1,0); //向量快速乘？【逃
    for(;p;p>>=1,x=x*x) if(p&1) s=s*x; return s.x;
}
inline int Sqrt(int x){ if(!x) return 0; // 0 的情况返回 0 就好了
    if(qpow(x,(mod-1)>>1)==mod-1) return -1; // 无解返回 -1
    while(1){ a=mul(rand(),rand()),w=dec(mul(a,a),x);
        if(qpow(w,(mod-1)>>1)==mod-1) return qpow(cp(a,1),(mod+1)>>1);
    }
}

ll _qpow(ll x,ll k)
{
    ll res = 1;
    while(k)
    {
        if(k & 1)
            res = res * x % mod;
        x = x * x % mod;
        k /= 2;
    }
    return res;
}
map<ll,ll>mp[2];
ll BSGS(ll q,ll v,int kind)
{
    mp[0].clear();
    mp[1].clear();
    ll T = sqrt(mod) + 1;
    ll res = v;
    //mp[0][res] = -1;
    for(int i = 1;i <= T;i++)
    {
        res = res * q % mod;
        if(!mp[i&1][res])
            mp[i&1][res] = i;
    }
    ll temp = _qpow(q,T);
    ll temp2;
    res = temp;
    for(int i = 1;i <= T;i++)
    {
        temp2 = mp[i * T & 1 ^ kind][res];
        if(temp2)
        {
            //if(temp2 == -1)
            //temp2 = 0;
            return i * T - temp2;
        }
        res = res * temp % mod;
    }
    return -1;
}
int main()
{
    srand(time(NULL));
    int t;
    ll c,two,five,q;
    scanf("%d",&t);
    while(t--)
    {
        scanf("%lld%lld",&c,&mod);
        five = Sqrt(5);
        //printf("%lld\n",five);
        two = _qpow(2,mod - 2);
        c = c * five % mod;
        q = (five + 1) * two % mod;
        
        ll ans = 2e9 + 1,res;
        ll b,v;
        
        b = Sqrt((c * c + 4) % mod);
        if(b >= 0)
        {
            v = (c + b) % mod * two % mod;
            res = BSGS(q,v,0);
            if(res != -1)
                ans = min(ans,res);
            
            v = ((c - b) % mod + mod) % mod;
            v = v * two % mod;
            res = BSGS(q,v,0);
            if(res != -1)
                ans = min(ans,res);
        }
        
        b = Sqrt((c * c + mod - 4) % mod);
        if(b >= 0)
        {
            v = (c + b) % mod * two % mod;
            res = BSGS(q,v,1);
            if(res != -1)
                ans = min(ans,res);
            
            v = ((c - b) % mod + mod) % mod;
            v = v * two % mod;
            res = BSGS(q,v,1);
            if(res != -1)
                ans = min(ans,res);
        }
        
        if(ans < 2e9 + 1)
            printf("%lld\n",ans);
        else
            printf("-1\n",ans);
    }
    return 0;
}
