#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(PRJ_ENABLE_MPI)
#include <mpi.h>
#endif

#include "prj.h"

#if PRJ_DYNAMIC_GR

typedef struct prj_wave_complex {
    double re;
    double im;
} prj_wave_complex;

typedef struct prj_wave_vertex {
    double x[3];
    double theta;
    double phi;
    double weight;
} prj_wave_vertex;

typedef struct prj_wave_grid {
    size_t nvert;
    size_t ntri;
    prj_wave_vertex *vert;
    size_t *tri;
} prj_wave_grid;

typedef struct prj_wave_hash_slot {
    int used;
    int64_t q[3];
    size_t index;
} prj_wave_hash_slot;

static void prj_wave_fail(const char *message)
{
    fprintf(stderr, "%s\n", message);
#if defined(PRJ_ENABLE_MPI)
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
    exit(EXIT_FAILURE);
}

static int prj_wave_is_root(const prj_mpi *mpi)
{
    return mpi == 0 || mpi->rank == 0;
}

static int prj_wave_local_block(const prj_mpi *mpi, const prj_block *block)
{
    return block != 0 && block->id >= 0 && block->active == 1 && block->z4c != 0 &&
        (mpi == 0 || block->rank == mpi->rank);
}

static size_t prj_wave_checked_mul(size_t a, size_t b, const char *what)
{
    if (a != 0 && b > SIZE_MAX / a) {
        prj_wave_fail(what);
    }
    return a * b;
}

static size_t prj_wave_mode_count(int lmax)
{
    size_t lp1;
    size_t square;

    if (lmax < 2) {
        prj_wave_fail("prj_z4c_extract_wave: z4c_wave_lmax must be at least 2");
    }
    lp1 = (size_t)lmax + 1U;
    square = prj_wave_checked_mul(lp1, lp1,
        "prj_z4c_extract_wave: z4c_wave_lmax overflows the mode count");
    return square - 4U;
}

static size_t prj_wave_mode_index(int l, int m)
{
    return (size_t)l * (size_t)l - 4U + (size_t)(m + l);
}

/* Physical ADM field sampled at a Z4c cell (including Z4c ghost cells). */
static double prj_wave_adm_field(const prj_mesh *mesh, const prj_block *block,
    int i, int j, int k, int field, int a, int b)
{
    double gamma[3][3];
    double K_dd[3][3];

    if (!prj_z4c_load_adm_cell(mesh, block, 0, i, j, k, gamma, K_dd)) {
        prj_wave_fail("prj_z4c_extract_wave: failed to reconstruct ADM geometry");
    }
    return field == 0 ? gamma[a][b] : K_dd[a][b];
}

static double prj_wave_fd1(const prj_mesh *mesh, const prj_block *block,
    int field, int a, int b, int dir, int i, int j, int k)
{
#if PRJ_NGHOST_Z4C == 2
    static const int off[4] = {-2, -1, 1, 2};
    static const double c[4] = {1.0/12.0, -2.0/3.0, 2.0/3.0, -1.0/12.0};
    const int n = 4;
#elif PRJ_NGHOST_Z4C == 4
    static const int off[8] = {-4, -3, -2, -1, 1, 2, 3, 4};
    static const double c[8] = {
        1.0/280.0, -4.0/105.0, 1.0/5.0, -4.0/5.0,
        4.0/5.0, -1.0/5.0, 4.0/105.0, -1.0/280.0
    };
    const int n = 8;
#else
#error "GW extraction requires PRJ_NGHOST_Z4C == 2 or 4"
#endif
    double sum = 0.0;
    int p;

    for (p = 0; p < n; ++p) {
        int ii = i + (dir == 0 ? off[p] : 0);
        int jj = j + (dir == 1 ? off[p] : 0);
        int kk = k + (dir == 2 ? off[p] : 0);
        sum += c[p] * prj_wave_adm_field(mesh, block, ii, jj, kk, field, a, b);
    }
    return sum / block->dx[dir];
}

static double prj_wave_fd2(const prj_mesh *mesh, const prj_block *block,
    int field, int a, int b, int dir, int i, int j, int k)
{
#if PRJ_NGHOST_Z4C == 2
    static const int off[5] = {-2, -1, 0, 1, 2};
    static const double c[5] = {-1.0/12.0, 4.0/3.0, -5.0/2.0, 4.0/3.0, -1.0/12.0};
    const int n = 5;
#elif PRJ_NGHOST_Z4C == 4
    static const int off[9] = {-4, -3, -2, -1, 0, 1, 2, 3, 4};
    static const double c[9] = {
        -1.0/560.0, 8.0/315.0, -1.0/5.0, 8.0/5.0, -205.0/72.0,
        8.0/5.0, -1.0/5.0, 8.0/315.0, -1.0/560.0
    };
    const int n = 9;
#endif
    double sum = 0.0;
    int p;

    for (p = 0; p < n; ++p) {
        int ii = i + (dir == 0 ? off[p] : 0);
        int jj = j + (dir == 1 ? off[p] : 0);
        int kk = k + (dir == 2 ? off[p] : 0);
        sum += c[p] * prj_wave_adm_field(mesh, block, ii, jj, kk, field, a, b);
    }
    return sum / (block->dx[dir] * block->dx[dir]);
}

static double prj_wave_fdxy(const prj_mesh *mesh, const prj_block *block,
    int field, int a, int b, int d1, int d2, int i, int j, int k)
{
#if PRJ_NGHOST_Z4C == 2
    static const int off[4] = {-2, -1, 1, 2};
    static const double c[4] = {1.0/12.0, -2.0/3.0, 2.0/3.0, -1.0/12.0};
    const int n = 4;
#elif PRJ_NGHOST_Z4C == 4
    static const int off[8] = {-4, -3, -2, -1, 1, 2, 3, 4};
    static const double c[8] = {
        1.0/280.0, -4.0/105.0, 1.0/5.0, -4.0/5.0,
        4.0/5.0, -1.0/5.0, 4.0/105.0, -1.0/280.0
    };
    const int n = 8;
#endif
    double sum = 0.0;
    int p, q;

    for (p = 0; p < n; ++p) {
        for (q = 0; q < n; ++q) {
            int ii = i + (d1 == 0 ? off[p] : 0) + (d2 == 0 ? off[q] : 0);
            int jj = j + (d1 == 1 ? off[p] : 0) + (d2 == 1 ? off[q] : 0);
            int kk = k + (d1 == 2 ? off[p] : 0) + (d2 == 2 ? off[q] : 0);
            sum += c[p] * c[q] *
                prj_wave_adm_field(mesh, block, ii, jj, kk, field, a, b);
        }
    }
    return sum / (block->dx[d1] * block->dx[d2]);
}

static double prj_wave_det3(const double g[3][3])
{
    return g[0][0]*(g[1][1]*g[2][2] - g[1][2]*g[2][1])
        - g[0][1]*(g[1][0]*g[2][2] - g[1][2]*g[2][0])
        + g[0][2]*(g[1][0]*g[2][1] - g[1][1]*g[2][0]);
}

static void prj_wave_inv3(const double g[3][3], double u[3][3])
{
    double det = prj_wave_det3(g);

    if (!isfinite(det) || det <= DBL_MIN) {
        prj_wave_fail("prj_z4c_extract_wave: non-positive physical metric determinant");
    }
    u[0][0] = (g[1][1]*g[2][2] - g[1][2]*g[2][1])/det;
    u[0][1] = (g[0][2]*g[2][1] - g[0][1]*g[2][2])/det;
    u[0][2] = (g[0][1]*g[1][2] - g[0][2]*g[1][1])/det;
    u[1][0] = u[0][1];
    u[1][1] = (g[0][0]*g[2][2] - g[0][2]*g[2][0])/det;
    u[1][2] = (g[0][2]*g[1][0] - g[0][0]*g[1][2])/det;
    u[2][0] = u[0][2];
    u[2][1] = u[1][2];
    u[2][2] = (g[0][0]*g[1][1] - g[0][1]*g[1][0])/det;
}

static double prj_wave_dot(const double g[3][3], const double a[3], const double b[3])
{
    double out = 0.0;
    int i, j;
    for (i = 0; i < 3; ++i) for (j = 0; j < 3; ++j) out += g[i][j]*a[i]*b[j];
    return out;
}

static void prj_wave_normalize(const double g[3][3], double v[3])
{
    double n = prj_wave_dot(g, v, v);
    int a;
    if (!isfinite(n) || n <= DBL_MIN) prj_wave_fail("prj_z4c_extract_wave: degenerate wave tetrad");
    n = sqrt(n);
    for (a = 0; a < 3; ++a) v[a] /= n;
}

static prj_wave_complex prj_wave_psi4_cell(const prj_mesh *mesh,
    const prj_block *block, int i, int j, int k)
{
    double g[3][3], gu[3][3], Kdd[3][3], Kud[3][3];
    double dg[3][3][3], dK[3][3][3], ddg[3][3][3][3];
    double G1[3][3][3], G2[3][3][3], Ric[3][3], DK[3][3][3];
    double R3[3][3][3][3], R4[3][3][3][3], Rn[3][3][3], Rnn[3][3];
    double rv[3], tv[3], pv[3];
    double R = 0.0, K = 0.0;
    prj_wave_complex out = {0.0, 0.0};
    int a,b,c,d,e;

    if (!prj_z4c_load_adm_cell(mesh, block, 0, i, j, k, g, Kdd))
        prj_wave_fail("prj_z4c_extract_wave: failed to load central ADM geometry");
    prj_wave_inv3((const double (*)[3])g, gu);
    memset(dg,0,sizeof(dg)); memset(dK,0,sizeof(dK)); memset(ddg,0,sizeof(ddg));
    memset(G1,0,sizeof(G1)); memset(G2,0,sizeof(G2)); memset(Ric,0,sizeof(Ric));
    memset(DK,0,sizeof(DK)); memset(R3,0,sizeof(R3)); memset(R4,0,sizeof(R4));
    memset(Rn,0,sizeof(Rn)); memset(Rnn,0,sizeof(Rnn)); memset(Kud,0,sizeof(Kud));

    for (c=0;c<3;++c) for (a=0;a<3;++a) for (b=0;b<3;++b) {
        dg[c][a][b]=prj_wave_fd1(mesh,block,0,a,b,c,i,j,k);
        dK[c][a][b]=prj_wave_fd1(mesh,block,1,a,b,c,i,j,k);
    }
    for (c=0;c<3;++c) for (d=c;d<3;++d) for (a=0;a<3;++a) for (b=0;b<3;++b) {
        double v=c==d?prj_wave_fd2(mesh,block,0,a,b,c,i,j,k):
            prj_wave_fdxy(mesh,block,0,a,b,c,d,i,j,k);
        ddg[c][d][a][b]=ddg[d][c][a][b]=v;
    }
    for (c=0;c<3;++c) for (a=0;a<3;++a) for (b=0;b<3;++b) {
        G1[c][a][b]=0.5*(dg[a][b][c]+dg[b][a][c]-dg[c][a][b]);
        for (d=0;d<3;++d) G2[c][a][b]+=gu[c][d]*G1[d][a][b];
    }
    for (a=0;a<3;++a) for (b=0;b<3;++b) {
        for (c=0;c<3;++c) for (d=0;d<3;++d) {
            Ric[a][b]+=0.5*gu[c][d]*(-ddg[c][d][a][b]-ddg[a][b][c][d]
                +ddg[a][c][b][d]+ddg[b][c][a][d]);
            for (e=0;e<3;++e) {
                Ric[a][b]+=gu[c][d]*(G2[e][a][c]*G1[e][b][d]
                    -G2[e][a][b]*G1[e][c][d]);
            }
        }
        R+=gu[a][b]*Ric[a][b];
        for (c=0;c<3;++c) Kud[a][b]+=gu[a][c]*Kdd[c][b];
    }
    for (a=0;a<3;++a) K+=Kud[a][a];
    for (c=0;c<3;++c) for (a=0;a<3;++a) for (b=0;b<3;++b) {
        DK[c][a][b]=dK[c][a][b];
        for (d=0;d<3;++d) DK[c][a][b]-=G2[d][c][a]*Kdd[d][b]+G2[d][c][b]*Kdd[a][d];
    }
    for (a=0;a<3;++a) for (b=0;b<3;++b) for (c=0;c<3;++c) for (d=0;d<3;++d) {
        R3[a][b][c][d]=g[a][c]*Ric[b][d]+g[b][d]*Ric[a][c]
            -g[a][d]*Ric[b][c]-g[b][c]*Ric[a][d]
            -0.5*R*g[a][c]*g[b][d]+0.5*R*g[a][d]*g[b][c];
        R4[a][b][c][d]=R3[a][b][c][d]+Kdd[a][c]*Kdd[b][d]-Kdd[a][d]*Kdd[b][c];
    }
    for (a=0;a<3;++a) for (b=0;b<3;++b) for (c=0;c<3;++c)
        Rn[a][b][c]=-(DK[c][a][b]-DK[b][a][c]);
    for (a=0;a<3;++a) for (b=0;b<3;++b) {
        Rnn[a][b]=Ric[a][b]+K*Kdd[a][b];
        for (c=0;c<3;++c) for (d=0;d<3;++d) Rnn[a][b]-=gu[c][d]*Kdd[a][c]*Kdd[d][b];
    }

    rv[0]=block->xmin[0]+((double)i+0.5)*block->dx[0];
    rv[1]=block->xmin[1]+((double)j+0.5)*block->dx[1];
    rv[2]=block->xmin[2]+((double)k+0.5)*block->dx[2];
    if (rv[0]*rv[0]+rv[1]*rv[1] < 1.0e-20) rv[0]+=1.0e-10*block->dx[0];
    pv[0]=-rv[1]; pv[1]=rv[0]; pv[2]=0.0;
    tv[0]=rv[0]*rv[2]; tv[1]=rv[1]*rv[2]; tv[2]=-(rv[0]*rv[0]+rv[1]*rv[1]);
    prj_wave_normalize((const double (*)[3])g,pv);
    { double q=prj_wave_dot((const double (*)[3])g,pv,rv); for(a=0;a<3;++a) rv[a]-=q*pv[a]; }
    prj_wave_normalize((const double (*)[3])g,rv);
    { double q1=prj_wave_dot((const double (*)[3])g,pv,tv);
      double q2=prj_wave_dot((const double (*)[3])g,rv,tv);
      for(a=0;a<3;++a) tv[a]-=q1*pv[a]+q2*rv[a]; }
    prj_wave_normalize((const double (*)[3])g,tv);

    for (a=0;a<3;++a) for (b=0;b<3;++b) {
        double qr=tv[a]*tv[b]-pv[a]*pv[b];
        double qi=-(tv[a]*pv[b]+pv[a]*tv[b]);
        out.re-=0.25*Rnn[a][b]*qr; out.im-=0.25*Rnn[a][b]*qi;
        for(c=0;c<3;++c) {
            out.re+=0.5*Rn[a][c][b]*rv[c]*qr; out.im+=0.5*Rn[a][c][b]*rv[c]*qi;
            for(d=0;d<3;++d) {
                out.re-=0.25*R4[d][a][c][b]*rv[d]*rv[c]*qr;
                out.im-=0.25*R4[d][a][c][b]*rv[d]*rv[c]*qi;
            }
        }
    }
    if (!isfinite(out.re)||!isfinite(out.im)) prj_wave_fail("prj_z4c_extract_wave: non-finite Psi4");
    return out;
}

static uint64_t prj_wave_hash_key(int64_t x, int64_t y, int64_t z)
{
    uint64_t h=1469598103934665603ULL;
    h=(h^(uint64_t)x)*1099511628211ULL; h=(h^(uint64_t)y)*1099511628211ULL;
    return (h^(uint64_t)z)*1099511628211ULL;
}

static size_t prj_wave_hash_insert(prj_wave_hash_slot *tab, size_t cap,
    prj_wave_vertex *vert, size_t *nvert, const double p[3])
{
    int64_t q[3]; size_t slot; int a;
    for(a=0;a<3;++a) q[a]=(int64_t)llround(p[a]*1.0e13);
    slot=(size_t)(prj_wave_hash_key(q[0],q[1],q[2])&(uint64_t)(cap-1U));
    while(tab[slot].used) {
        if(tab[slot].q[0]==q[0]&&tab[slot].q[1]==q[1]&&tab[slot].q[2]==q[2]) return tab[slot].index;
        slot=(slot+1U)&(cap-1U);
    }
    tab[slot].used=1; memcpy(tab[slot].q,q,sizeof(q)); tab[slot].index=*nvert;
    memcpy(vert[*nvert].x,p,3U*sizeof(double)); vert[*nvert].weight=0.0;
    *nvert+=1U; return tab[slot].index;
}

static double prj_wave_euclid_dot(const double a[3], const double b[3])
{ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

static void prj_wave_cross(const double a[3], const double b[3], double c[3])
{ c[0]=a[1]*b[2]-a[2]*b[1]; c[1]=a[2]*b[0]-a[0]*b[2]; c[2]=a[0]*b[1]-a[1]*b[0]; }

static void prj_wave_unit(double p[3])
{ double n=sqrt(prj_wave_euclid_dot(p,p)); p[0]/=n;p[1]/=n;p[2]/=n; }

static double prj_wave_solid_triangle(const double a[3],const double b[3],const double c[3])
{
    double bx[3]; double det;
    prj_wave_cross(b,c,bx); det=fabs(prj_wave_euclid_dot(a,bx));
    return 2.0*atan2(det,1.0+prj_wave_euclid_dot(a,b)+prj_wave_euclid_dot(b,c)+prj_wave_euclid_dot(c,a));
}

static size_t prj_wave_tri_local(int n,int i,int j)
{ return (size_t)i*(size_t)(n+1)-(size_t)i*(size_t)(i-1)/2U+(size_t)j; }

static int prj_wave_angle_cmp(const void *aa,const void *bb)
{
    const double *a=(const double*)aa,*b=(const double*)bb;
    return a[0]<b[0]?-1:(a[0]>b[0]?1:0);
}

static prj_wave_grid prj_wave_make_grid(int nq,double radius)
{
    static const int faces[20][3]={{0,11,5},{0,5,1},{0,1,7},{0,7,10},{0,10,11},
      {1,5,9},{5,11,4},{11,10,2},{10,7,6},{7,1,8},{3,9,4},{3,4,2},{3,2,6},
      {3,6,8},{3,8,9},{4,9,5},{2,4,11},{6,2,10},{8,6,7},{9,8,1}};
    double base[12][3]; double phi=(1.0+sqrt(5.0))/2.0;
    prj_wave_grid grid={0,0,0,0}; size_t expected,ntri,cap=1,local_count;
    prj_wave_hash_slot *hash; size_t f,tpos=0; int i,j;
    const double raw[12][3]={{-1,phi,0},{1,phi,0},{-1,-phi,0},{1,-phi,0},
      {0,-1,phi},{0,1,phi},{0,-1,-phi},{0,1,-phi},{phi,0,-1},{phi,0,1},
      {-phi,0,-1},{-phi,0,1}};
    if(nq<1) prj_wave_fail("prj_z4c_extract_wave: invalid geodesic level");
    if((size_t)nq>sqrt((double)(SIZE_MAX/20U))) prj_wave_fail("prj_z4c_extract_wave: geodesic level overflow");
    expected=10U*(size_t)nq*(size_t)nq+2U; ntri=20U*(size_t)nq*(size_t)nq;
    while(cap<2U*expected) { if(cap>SIZE_MAX/2U) prj_wave_fail("prj_z4c_extract_wave: geodesic hash overflow"); cap*=2U; }
    grid.vert=(prj_wave_vertex*)prj_calloc(expected,sizeof(*grid.vert));
    grid.tri=(size_t*)prj_calloc(prj_wave_checked_mul(3U,ntri,"prj_z4c_extract_wave: triangle overflow"),sizeof(*grid.tri));
    hash=(prj_wave_hash_slot*)prj_calloc(cap,sizeof(*hash));
    local_count=prj_wave_checked_mul((size_t)nq+1U,(size_t)nq+2U,
        "prj_z4c_extract_wave: local geodesic size overflow")/2U;
    if(!grid.vert||!grid.tri||!hash) prj_wave_fail("prj_z4c_extract_wave: geodesic allocation failed");
    memcpy(base,raw,sizeof(base)); for(i=0;i<12;++i) prj_wave_unit(base[i]);
    for(f=0;f<20U;++f) {
        size_t *lid=(size_t*)prj_calloc(local_count,sizeof(*lid));
        if(!lid) prj_wave_fail("prj_z4c_extract_wave: local geodesic allocation failed");
        for(i=0;i<=nq;++i) for(j=0;j<=nq-i;++j) {
            double p[3]; int k=nq-i-j,a;
            for(a=0;a<3;++a) p[a]=((double)k*base[faces[f][0]][a]+(double)i*base[faces[f][1]][a]+(double)j*base[faces[f][2]][a])/(double)nq;
            prj_wave_unit(p); lid[prj_wave_tri_local(nq,i,j)]=prj_wave_hash_insert(hash,cap,grid.vert,&grid.nvert,p);
        }
        for(i=0;i<nq;++i) for(j=0;j<nq-i;++j) {
            grid.tri[tpos++]=lid[prj_wave_tri_local(nq,i,j)];
            grid.tri[tpos++]=lid[prj_wave_tri_local(nq,i+1,j)];
            grid.tri[tpos++]=lid[prj_wave_tri_local(nq,i,j+1)];
            if(j<nq-i-1) {
                grid.tri[tpos++]=lid[prj_wave_tri_local(nq,i+1,j)];
                grid.tri[tpos++]=lid[prj_wave_tri_local(nq,i+1,j+1)];
                grid.tri[tpos++]=lid[prj_wave_tri_local(nq,i,j+1)];
            }
        }
        free(lid);
    }
    free(hash); grid.ntri=tpos/3U;
    if(grid.nvert!=expected||grid.ntri!=ntri) prj_wave_fail("prj_z4c_extract_wave: invalid geodesic topology");
    {
        size_t *degree=(size_t*)prj_calloc(grid.nvert,sizeof(*degree));
        size_t *offset,*cursor; double (*centers)[3]; size_t *incident;
        if(!degree) prj_wave_fail("prj_z4c_extract_wave: degree allocation failed");
        for(f=0;f<grid.ntri;++f) for(i=0;i<3;++i) degree[grid.tri[3U*f+(size_t)i]]++;
        offset=(size_t*)prj_calloc(grid.nvert+1U,sizeof(*offset)); cursor=(size_t*)prj_calloc(grid.nvert,sizeof(*cursor));
        centers=(double(*)[3])prj_calloc(grid.ntri,sizeof(*centers));
        incident=(size_t*)prj_calloc(prj_wave_checked_mul(3U,grid.ntri,
            "prj_z4c_extract_wave: incidence size overflow"),sizeof(*incident));
        if(!offset||!cursor||!centers||!incident) prj_wave_fail("prj_z4c_extract_wave: dual-grid allocation failed");
        for(f=0;f<grid.nvert;++f) offset[f+1U]=offset[f]+degree[f]; memcpy(cursor,offset,grid.nvert*sizeof(*cursor));
        for(f=0;f<grid.ntri;++f) {
            const double *a=grid.vert[grid.tri[3U*f]].x,*b=grid.vert[grid.tri[3U*f+1U]].x,*c=grid.vert[grid.tri[3U*f+2U]].x;
            double ba[3]={b[0]-a[0],b[1]-a[1],b[2]-a[2]},ca[3]={c[0]-a[0],c[1]-a[1],c[2]-a[2]};
            prj_wave_cross(ba,ca,centers[f]); prj_wave_unit(centers[f]);
            if(prj_wave_euclid_dot(centers[f],a)<0) for(i=0;i<3;++i) centers[f][i]*=-1.0;
            for(i=0;i<3;++i) { size_t v=grid.tri[3U*f+(size_t)i]; incident[cursor[v]++]=f; }
        }
        for(f=0;f<grid.nvert;++f) {
            double ref[3]={0,0,1},e1[3],e2[3],items[12][2]; size_t q,n=degree[f];
            if(fabs(grid.vert[f].x[2])>0.9) { ref[0]=1;ref[2]=0; }
            prj_wave_cross(ref,grid.vert[f].x,e1);prj_wave_unit(e1);prj_wave_cross(grid.vert[f].x,e1,e2);
            if(n>12U) prj_wave_fail("prj_z4c_extract_wave: unexpected geodesic degree");
            for(q=0;q<n;++q) { size_t ti=incident[offset[f]+q]; items[q][0]=atan2(prj_wave_euclid_dot(centers[ti],e2),prj_wave_euclid_dot(centers[ti],e1)); items[q][1]=(double)ti; }
            qsort(items,n,sizeof(items[0]),prj_wave_angle_cmp);
            for(q=0;q<n;++q) { size_t t0=(size_t)items[q][1],t1=(size_t)items[(q+1U)%n][1]; grid.vert[f].weight+=prj_wave_solid_triangle(grid.vert[f].x,centers[t0],centers[t1]); }
        }
        free(degree);free(offset);free(cursor);free(centers);free(incident);
    }
    { double sum=0; for(f=0;f<grid.nvert;++f) { double *p=grid.vert[f].x; grid.vert[f].theta=acos(fmax(-1.0,fmin(1.0,p[2])));grid.vert[f].phi=atan2(p[1],p[0]);if(grid.vert[f].phi<0)grid.vert[f].phi+=2*M_PI;sum+=grid.vert[f].weight;p[0]*=radius;p[1]*=radius;p[2]*=radius; }
      if(fabs(sum-4.0*M_PI)>1.0e-10*fmax(1.0,4.0*M_PI)) prj_wave_fail("prj_z4c_extract_wave: geodesic weights do not sum to 4 pi"); }
    return grid;
}

static void prj_wave_free_grid(prj_wave_grid *grid)
{ if(grid){free(grid->vert);free(grid->tri);memset(grid,0,sizeof(*grid));} }

static int prj_wave_box_intersects_sphere(const prj_block *b,double r)
{
    double dmin=0,dmax=0;int a;
    for(a=0;a<3;++a){double lo=b->xmin[a],hi=b->xmax[a],mn=0,mx=fmax(fabs(lo),fabs(hi));if(0<lo)mn=lo;else if(0>hi)mn=-hi;dmin+=mn*mn;dmax+=mx*mx;}
    return r*r>=dmin*(1.0-1e-13)&&r*r<=dmax*(1.0+1e-13);
}

static int prj_wave_choose_nq(const prj_mesh *mesh,double radius,int lmax)
{
    double dx=HUGE_VAL,base,need;int bidx;size_t modes=prj_wave_mode_count(lmax);
    for(bidx=0;bidx<mesh->nblocks;++bidx){const prj_block*b=&mesh->blocks[bidx];if(b->id>=0&&b->active==1&&prj_wave_box_intersects_sphere(b,radius)){double d=fmin(b->dx[0],fmin(b->dx[1],b->dx[2]));if(d<dx)dx=d;}}
    if(!isfinite(dx)||dx<=0)prj_wave_fail("prj_z4c_extract_wave: extraction sphere does not intersect the active mesh");
    base=M_PI*radius*radius/(dx*dx)-2.0;if(base<1)base=1;need=((double)modes-2.0)/10.0;if(need>base)base=need;
    base=ceil(sqrt(base));if(!isfinite(base)||base>(double)INT_MAX)prj_wave_fail("prj_z4c_extract_wave: geodesic level exceeds integer range");return (int)base;
}

static int prj_wave_block_contains(const prj_mesh *mesh,const prj_block*b,const double x[3])
{
    const double domhi[3]={mesh->coord.x1max,mesh->coord.x2max,mesh->coord.x3max};int a;
    for(a=0;a<3;++a){double tol=1e-12*fmax(1.0,fabs(b->dx[a]));if(x[a]<b->xmin[a]-tol)return 0;if(x[a]>=b->xmax[a]-tol&&fabs(b->xmax[a]-domhi[a])>tol)return 0;if(x[a]>b->xmax[a]+tol)return 0;}return 1;
}

static void prj_wave_interp_weights(const prj_block*b,int dir,double x,int n,int*base,double*w)
{
    double t=(x-(b->xmin[dir]+0.5*b->dx[dir]))/b->dx[dir];int q,p,ib=(int)floor(t)-(n/2-1);
    if(ib<0)ib=0;if(ib>PRJ_BLOCK_SIZE-n)ib=PRJ_BLOCK_SIZE-n;*base=ib;
    for(q=0;q<n;++q){double v=1;for(p=0;p<n;++p)if(p!=q)v*=(t-(double)(ib+p))/(double)(q-p);w[q]=v;}
}

static prj_wave_complex prj_wave_interp_block(const prj_block*b,const prj_wave_complex*data,const double x[3])
{
    enum{N=2*PRJ_NGHOST_Z4C};double wx[N],wy[N],wz[N];int bi,bj,bk,i,j,k;prj_wave_complex o={0,0};
    prj_wave_interp_weights(b,0,x[0],N,&bi,wx);prj_wave_interp_weights(b,1,x[1],N,&bj,wy);prj_wave_interp_weights(b,2,x[2],N,&bk,wz);
    for(i=0;i<N;++i)for(j=0;j<N;++j)for(k=0;k<N;++k){double w=wx[i]*wy[j]*wz[k];const prj_wave_complex*v=&data[LIDX(bi+i,bj+j,bk+k)];o.re+=w*v->re;o.im+=w*v->im;}return o;
}

static double prj_wave_log_power(double x,int p)
{ if(p==0)return 0.0;if(x<=0)return -HUGE_VAL;return (double)p*log(x); }

static void prj_wave_swsh(int l,int m,double theta,double phi,double *yr,double *yi)
{
    int k,k1=m-2>0?m-2:0,k2=l+m<l-2?l+m:l-2;double ct=cos(0.5*theta),st=sin(0.5*theta),maxlog=-HUGE_VAL,sum=0;
    double norm=0.5*(lgamma((double)(l+m+1))+lgamma((double)(l-m+1))+lgamma((double)(l-1))+lgamma((double)(l+3)));
    for(k=k1;k<=k2;++k){int pc=2*l+m-2-2*k,ps=2*k-m+2;double lt=norm-lgamma((double)(l+m-k+1))-lgamma((double)(l-k-1))-lgamma((double)(k+1))-lgamma((double)(k-m+3))+prj_wave_log_power(ct,pc)+prj_wave_log_power(st,ps);if(lt>maxlog)maxlog=lt;}
    if(!isfinite(maxlog)){*yr=*yi=0;return;}
    for(k=k1;k<=k2;++k){int pc=2*l+m-2-2*k,ps=2*k-m+2;double lt=norm-lgamma((double)(l+m-k+1))-lgamma((double)(l-k-1))-lgamma((double)(k+1))-lgamma((double)(k-m+3))+prj_wave_log_power(ct,pc)+prj_wave_log_power(st,ps);sum+=(k&1?-1.0:1.0)*exp(lt-maxlog);}
    sum*=exp(maxlog)*sqrt((2.0*(double)l+1.0)/(4.0*M_PI));*yr=sum*cos((double)m*phi);*yi=sum*sin((double)m*phi);
}

static void prj_wave_reduce_modes(const prj_mpi*mpi,double*local,double*global,size_t n)
{
#if defined(PRJ_ENABLE_MPI)
    if(mpi&&mpi->totrank>1){size_t p=0;while(p<n){int chunk=(n-p)>(size_t)INT_MAX?INT_MAX:(int)(n-p);MPI_Reduce(local+p,global+p,chunk,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);p+=(size_t)chunk;}return;}
#else
    (void)mpi;
#endif
    memcpy(global,local,n*sizeof(*global));
}

static void prj_wave_reduce_owners(const prj_mpi *mpi, int *local, int *global,
    size_t n)
{
#if defined(PRJ_ENABLE_MPI)
    if (mpi && mpi->totrank > 1) {
        size_t p = 0;
        while (p < n) {
            int chunk = n - p > (size_t)INT_MAX ? INT_MAX : (int)(n - p);
            MPI_Allreduce(local + p, global + p, chunk, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
            p += (size_t)chunk;
        }
        return;
    }
#else
    (void)mpi;
#endif
    memcpy(global, local, n * sizeof(*global));
}

void prj_z4c_wave_init(const prj_mpi *mpi)
{
    FILE*fp;if(!prj_wave_is_root(mpi))return;fp=fopen("output/z4c_wave.txt","w");if(!fp){fprintf(stderr,"failed to truncate output/z4c_wave.txt: %s\n",strerror(errno));return;}fprintf(fp,"# time_s radius_cm l m psi4_real_cm^-2 psi4_imag_cm^-2\n");fclose(fp);
}

void prj_z4c_extract_wave(const prj_mesh *mesh,const prj_mpi *mpi,double time_seconds)
{
    int lmax,nq,bidx,l,m;double radius;size_t nmodes,nvals,v;prj_wave_grid grid;double *local,*global;int *owners,*owners_global;
    if(!mesh)prj_wave_fail("prj_z4c_extract_wave: mesh is null");radius=mesh->z4c_params.wave_radii_cm;lmax=mesh->z4c_params.wave_lmax;
    if(!isfinite(radius)||radius<=0)prj_wave_fail("prj_z4c_extract_wave: invalid extraction radius");nmodes=prj_wave_mode_count(lmax);nvals=prj_wave_checked_mul(2U,nmodes,"prj_z4c_extract_wave: mode buffer overflow");nq=prj_wave_choose_nq(mesh,radius,lmax);grid=prj_wave_make_grid(nq,radius);
    local=(double*)prj_calloc(nvals,sizeof(*local));global=(double*)prj_calloc(nvals,sizeof(*global));owners=(int*)prj_calloc(grid.nvert,sizeof(*owners));owners_global=(int*)prj_calloc(grid.nvert,sizeof(*owners_global));
    if(!local||!global||!owners||!owners_global)prj_wave_fail("prj_z4c_extract_wave: output allocation failed");
    for(bidx=0;bidx<mesh->nblocks;++bidx){const prj_block*b=&mesh->blocks[bidx];prj_wave_complex*cell;int i,j,k;if(!prj_wave_local_block(mpi,b))continue;cell=(prj_wave_complex*)prj_calloc(PRJ_BLOCK_NCELLS,sizeof(*cell));if(!cell)prj_wave_fail("prj_z4c_extract_wave: cell allocation failed");for(i=0;i<PRJ_BLOCK_SIZE;++i)for(j=0;j<PRJ_BLOCK_SIZE;++j)for(k=0;k<PRJ_BLOCK_SIZE;++k)cell[LIDX(i,j,k)]=prj_wave_psi4_cell(mesh,b,i,j,k);
        for(v=0;v<grid.nvert;++v)if(prj_wave_block_contains(mesh,b,grid.vert[v].x)){prj_wave_complex p=prj_wave_interp_block(b,cell,grid.vert[v].x);owners[v]++;for(l=2;l<=lmax;++l)for(m=-l;m<=l;++m){double yr,yi,w=grid.vert[v].weight;size_t q=prj_wave_mode_index(l,m);prj_wave_swsh(l,m,grid.vert[v].theta,grid.vert[v].phi,&yr,&yi);local[2U*q]+=w*(p.re*yr+p.im*yi);local[2U*q+1U]+=w*(p.im*yr-p.re*yi);}}free(cell);}
    prj_wave_reduce_owners(mpi, owners, owners_global, grid.nvert);
    for(v=0;v<grid.nvert;++v)if(owners_global[v]!=1)prj_wave_fail("prj_z4c_extract_wave: extraction vertex does not have exactly one owner");
    prj_wave_reduce_modes(mpi,local,global,nvals);
    if(prj_wave_is_root(mpi)){FILE*fp=fopen("output/z4c_wave.txt","a");if(!fp)prj_wave_fail("prj_z4c_extract_wave: failed to open output/z4c_wave.txt");for(l=2;l<=lmax;++l)for(m=-l;m<=l;++m){size_t q=prj_wave_mode_index(l,m);fprintf(fp,"%.17e %.17e %d %d %.17e %.17e\n",time_seconds,radius,l,m,global[2U*q],global[2U*q+1U]);}fclose(fp);}
    free(local);free(global);free(owners);free(owners_global);prj_wave_free_grid(&grid);
}

int prj_z4c_wave_self_test(void)
{
    prj_wave_grid grid = prj_wave_make_grid(12, 1.0);
    double norm = 0.0;
    double cross_re = 0.0;
    double yr, yi;
    size_t v;
    int status = 0;

    if (grid.nvert != 1442U || grid.ntri != 2880U) {
        status = 1;
    }
    for (v = 0; v < grid.nvert; ++v) {
        double zr, zi;
        if (!isfinite(grid.vert[v].weight) || grid.vert[v].weight <= 0.0) {
            status = 2;
            break;
        }
        prj_wave_swsh(2, 2, grid.vert[v].theta, grid.vert[v].phi, &yr, &yi);
        prj_wave_swsh(2, 1, grid.vert[v].theta, grid.vert[v].phi, &zr, &zi);
        norm += grid.vert[v].weight * (yr * yr + yi * yi);
        cross_re += grid.vert[v].weight * (yr * zr + yi * zi);
    }
    if (fabs(norm - 1.0) > 5.0e-4 || fabs(cross_re) > 5.0e-4) {
        status = 3;
    }
    prj_wave_swsh(33, 17, 1.1, 0.7, &yr, &yi);
    if (!isfinite(yr) || !isfinite(yi) || prj_wave_mode_count(33) != 1152U) {
        status = 4;
    }
    prj_wave_free_grid(&grid);
    return status;
}

#else

void prj_z4c_wave_init(const prj_mpi *mpi) { (void)mpi; }
void prj_z4c_extract_wave(const prj_mesh *mesh,const prj_mpi *mpi,double time_seconds)
{ (void)mesh;(void)mpi;(void)time_seconds; }
int prj_z4c_wave_self_test(void) { return 0; }

#endif
