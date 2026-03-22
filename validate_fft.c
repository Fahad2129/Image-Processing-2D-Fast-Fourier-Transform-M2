/*
 * validate_fft.c  –  Reference validator for Milestone 2
 *
 * Compile natively (on your laptop):
 *   gcc -O0 -o validate validate_fft.c -lm && ./validate
 *
 * This prints the EXPECTED FFT output for each test size so you can
 * compare against the RISC-V assembly output.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846f

/* ── Taylor sin/cos (same as assembly version) ─────────────────────── */
static float wrap_angle(float x) {
    while (x >  (float)PI) x -= 2.0f*(float)PI;
    while (x < -(float)PI) x += 2.0f*(float)PI;
    return x;
}
static float my_sin(float x) {
    x = wrap_angle(x);
    float x2 = x*x;
    return x*(1.0f + x2*(-1.0f/6.0f
           + x2*( 1.0f/120.0f
           + x2*(-1.0f/5040.0f
           + x2*( 1.0f/362880.0f
           + x2*(-1.0f/39916800.0f))))));
}
static float my_cos(float x){ return my_sin(x + (float)PI/2.0f); }

/* ── Twiddle generation ─────────────────────────────────────────────── */
static void gen_twiddle(float *wr, float *wi, int n) {
    for (int k = 0; k < n/2; k++) {
        float a = -2.0f*(float)PI*k/n;
        wr[k] = my_cos(a);
        wi[k] = my_sin(a);
    }
}

/* ── Integer log2 ───────────────────────────────────────────────────── */
static int log2i(int n){ int l=0; while(n>1){n>>=1;l++;} return l; }

/* ── Bit reversal ───────────────────────────────────────────────────── */
static unsigned rev_bits(unsigned x, int logn) {
    unsigned r=0;
    for(int i=0;i<logn;i++){r=(r<<1)|(x&1);x>>=1;}
    return r;
}
static void bit_rev_arr(float *xr, float *xi, int n){
    int ln=log2i(n);
    for(int i=0;i<n;i++){
        unsigned j=rev_bits(i,ln);
        if(j>(unsigned)i){
            float t=xr[i];xr[i]=xr[j];xr[j]=t;
            t=xi[i];xi[i]=xi[j];xi[j]=t;
        }
    }
}

/* ── Iterative butterfly ────────────────────────────────────────────── */
static void butterfly_iter(float *xr, float *xi,
                            float *wr, float *wi, int n){
    int logn=log2i(n);
    for(int s=1;s<=logn;s++){
        int m=1<<s, half=m>>1, step=n/m;
        for(int k=0;k<n;k+=m)
            for(int j=0;j<half;j++){
                int ti=j*step;
                float Wr=wr[ti], Wi=wi[ti];
                float ur=xr[k+j], ui=xi[k+j];
                float vr=xr[k+j+half], vi=xi[k+j+half];
                float tr=Wr*vr-Wi*vi, t2=Wr*vi+Wi*vr;
                xr[k+j]=ur+tr; xi[k+j]=ui+t2;
                xr[k+j+half]=ur-tr; xi[k+j+half]=ui-t2;
            }
    }
}

/* ── Recursive butterfly helper ─────────────────────────────────────── */
static void rec_helper(float *xr, float *xi,
                        float *wr, float *wi, int n, int full_n){
    if(n<=1)return;
    int half=n/2, step=full_n/n;
    rec_helper(xr,xi,wr,wi,half,full_n);
    rec_helper(xr+half,xi+half,wr,wi,half,full_n);
    for(int k=0;k<half;k++){
        int ti=k*step;
        float Wr=wr[ti], Wi=wi[ti];
        float ur=xr[k], ui=xi[k];
        float vr=xr[k+half], vi=xi[k+half];
        float tr=Wr*vr-Wi*vi, t2=Wr*vi+Wi*vr;
        xr[k]=ur+tr; xi[k]=ui+t2;
        xr[k+half]=ur-tr; xi[k+half]=ui-t2;
    }
}

/* ── Public wrappers ─────────────────────────────────────────────────── */
static void fft_iter(float *xr, float *xi, float *wr, float *wi, int n){
    bit_rev_arr(xr,xi,n);
    butterfly_iter(xr,xi,wr,wi,n);
}
static void fft_rec(float *xr, float *xi, float *wr, float *wi, int n){
    bit_rev_arr(xr,xi,n);
    rec_helper(xr,xi,wr,wi,n,n);
}

/* ── Helpers ─────────────────────────────────────────────────────────── */
static void print_result(const char *label, float *xr, float *xi, int n){
    printf("\n%s (N=%d):\n", label, n);
    for(int i=0;i<n;i++)
        printf("  [%2d]  %9.4f + %9.4fi\n", i, xr[i], xi[i]);
}

#define N_MAX 32
static float xr[N_MAX], xi[N_MAX], wr[N_MAX], wi[N_MAX];

static void run(int n, const char *tag) {
    /* impulse: x[0]=1, rest 0 */
    memset(xr,0,sizeof(xr)); memset(xi,0,sizeof(xi));
    xr[0]=1.0f;
    gen_twiddle(wr,wi,n);
    float wr2[N_MAX], wi2[N_MAX];
    float xr2[N_MAX], xi2[N_MAX];
    memcpy(xr2,xr,n*sizeof(float));
    memcpy(xi2,xi,n*sizeof(float));
    memcpy(wr2,wr,(n/2)*sizeof(float));
    memcpy(wi2,wi,(n/2)*sizeof(float));

    fft_iter(xr,xi,wr,wi,n);
    char buf[64]; snprintf(buf,64,"Iterative %s",tag);
    print_result(buf,xr,xi,n);

    fft_rec(xr2,xi2,wr2,wi2,n);
    snprintf(buf,64,"Recursive %s",tag);
    print_result(buf,xr2,xi2,n);
}

int main(void){
    printf("=== FFT Reference Validator (C) ===\n");
    printf("Expected output: impulse input → all 1.0+0.0i\n");

    run(8,  "N=8");
    run(16, "N=16");
    run(32, "N=32");

    printf("\n=== Bit Reversal (N=8, log2=3) ===\n");
    for(int i=0;i<8;i++)
        printf("  bit_reverse(%d, 3) = %u\n", i, rev_bits(i,3));

    printf("\n=== Twiddle Factors (N=8) ===\n");
    gen_twiddle(wr,wi,8);
    for(int k=0;k<4;k++)
        printf("  W[%d] = %9.5f + %9.5fi\n", k, wr[k], wi[k]);

    return 0;
}
