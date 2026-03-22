/*
 * main.c  -  Driver for RISC-V assembly 1D FFT (Milestone 2)
 *
 * Build:
 *   riscv64-linux-gnu-gcc -static -O0 -g -o fft1d main.c fft_asm.s -lm
 * Run:
 *   qemu-riscv64 ./fft1d
 *
 */
#include <stdio.h>
#include <stdlib.h>   /* malloc, free */
#include <string.h>
#include <math.h>

/* ── Assembly function declarations ─────────────────────────────────── */
extern float        sin_approx(float x);
extern float        cos_approx(float x);
extern int          log2_int(int n);
extern unsigned int reverse_bits(unsigned int x, int logn);
extern void         generate_twiddle_factors(float *tr, float *ti, int n);
extern void         bit_reverse_array(float *xr, float *xi, int n);
extern void         butterfly_iterative(float *xr, float *xi,
                                        float *tr, float *ti, int n);
extern void         butterfly_recursive_helper(float *xr, float *xi,
                                               float *tr, float *ti,
                                               int n, int full_n);

/* ── FFT wrappers (call assembly) ────────────────────────────────────── */
void fft_1d_iterative(float *xr, float *xi,
                      float *tr, float *ti, int n) {
    bit_reverse_array(xr, xi, n);
    butterfly_iterative(xr, xi, tr, ti, n);
}

void fft_1d_recursive(float *xr, float *xi,
                      float *tr, float *ti, int n) {
    bit_reverse_array(xr, xi, n);
    butterfly_recursive_helper(xr, xi, tr, ti, n, n);
}

/* ── Cycle counter ───────────────────────────────────────────────────── */
static inline unsigned long long read_cycles(void) {
    unsigned long long c;
    __asm__ volatile ("rdcycle %0" : "=r"(c));
    return c;
}

/* ── Stack memory from assembly frame sizes ──────────────────────────
 *
 * ITERATIVE  butterfly_iterative frame (fft_asm.s):
 *   addi sp, sp, -80   →  80 bytes
 *   Saves: ra,s0-s8    →  10 regs × 8 bytes = 80 bytes  ✓
 *   bit_reverse_array is called inside (uses its own frame of 48 bytes)
 *   but the iterative loop itself does NOT recurse, so peak stack =
 *   butterfly_iterative(80) + bit_reverse_array(48) + log2_int(leaf,~16)
 *   = ~144 bytes.  We report 176 to include log2_int called from
 *   bit_reverse_array (adds one more small frame).  Fixed constant
 *   because there is NO recursion — depth is always 3 frames.
 *
 * RECURSIVE  butterfly_recursive_helper frame (fft_asm.s):
 *   addi sp, sp, -64   →  64 bytes per recursive call
 *   Saves: ra,s0-s6    →  8 regs × 8 bytes = 64 bytes  ✓
 *   Recursion depth = log2(N).
 *   Base overhead (bit_reverse_array + log2_int) ≈ 96 bytes (fixed).
 *   Total = 96 + 64 * log2(N)
 * ─────────────────────────────────────────────────────────────────── */
static int iter_stack(int n) {
    (void)n;
    /* Fixed depth: butterfly_iterative(80) + bit_reverse_array(48)
     * + log2_int inside bit_reverse_array(16) + log2_int in main frame(16)
     * + small call overhead padding  = 176 bytes (no recursion) */
    return 176;
}
static int recu_stack(int n) {
    int levels = 0, tmp = n;
    while (tmp > 1) { tmp >>= 1; levels++; }
    /* 96 = bit_reverse_array(48) + log2_int(16) + top-level overhead(32)
     * 64 = one recursive butterfly_recursive_helper frame              */
    return 96 + 64 * levels;
}

/* ── Print helpers ───────────────────────────────────────────────────── */
static void print_array(const char *label, float *xr, float *xi, int n) {
    printf("%s\nOutput:\n", label);
    for (int i = 0; i < n; i++)
        printf("  [%2d]  %8.4f + %8.4fi\n", i,
               (double)xr[i], (double)xi[i]);
}

/* ── Run one FFT test ────────────────────────────────────────────────── */
/*
 * tr / ti are HEAP-allocated by the caller (see main).
 * They must be at least n/2 floats in size.
 */
static unsigned long long run_test(const char *label,
                                   float *in_r, float *in_i,
                                   float *tr, float *ti,
                                   int n, int is_iter) {
    float xr[32] = {0}, xi[32] = {0};
    memcpy(xr, in_r, n * sizeof(float));
    memcpy(xi, in_i, n * sizeof(float));

    /* generate_twiddle_factors writes into heap-allocated tr/ti */
    generate_twiddle_factors(tr, ti, n);   /* ASSEMBLY */

    unsigned long long t0 = read_cycles();
    if (is_iter) fft_1d_iterative(xr, xi, tr, ti, n);
    else         fft_1d_recursive(xr, xi, tr, ti, n);
    unsigned long long t1 = read_cycles();

    print_array(label, xr, xi, n);
    return t1 - t0;
}

/* ── main ────────────────────────────────────────────────────────────── */
int main(void) {

    /* 1. log2_int demo */
    printf("\n=== LOG2_INT DEMO ===\n");
    printf("  log2_int(8)  = %d\n", log2_int(8));
    printf("  log2_int(16) = %d\n", log2_int(16));
    printf("  log2_int(32) = %d\n", log2_int(32));

    /* 2. Bit-reversal demo */
    printf("\n=== BIT-REVERSAL TEST (N=8, log2=3) ===\n");
    for (int i = 0; i < 8; i++)
        printf("  bit_reverse(%d, 3) = %u\n", i, reverse_bits(i, 3));

    /* 3. Twiddle factor demo
     * ----------------------------------------------------------------
     * HEAP ALLOCATION: twiddle arrays are now malloc'd on the heap
     * as required by the spec (Section 3.2).
     * For N=8 we need N/2 = 4 complex twiddle factors.
     * ---------------------------------------------------------------- */
    int demo_n = 8;
    float *tw_r = (float *)malloc((demo_n / 2) * sizeof(float));
    float *tw_i = (float *)malloc((demo_n / 2) * sizeof(float));
    if (!tw_r || !tw_i) { fprintf(stderr, "malloc failed\n"); return 1; }

    generate_twiddle_factors(tw_r, tw_i, demo_n);   /* ASSEMBLY */
    printf("\n=== TWIDDLE FACTORS (N=8) — computed by assembly ===\n");
    printf("  (Reference: W[k] = cos(-2*PI*k/8) + j*sin(-2*PI*k/8))\n");
    for (int k = 0; k < demo_n / 2; k++) {
        float a = -2.0f * 3.14159265358979f * k / (float)demo_n;
        float ref_r = (float)cos(a);
        float ref_i = (float)sin(a);
        printf("  W[%d] asm=(%8.5f + %8.5fi)  ref=(%8.5f + %8.5fi)\n",
               k, (double)tw_r[k], (double)tw_i[k],
               (double)ref_r, (double)ref_i);
    }
    free(tw_r);
    free(tw_i);

    /* 4. FFT tests
     * ----------------------------------------------------------------
     * HEAP ALLOCATION: one pair of twiddle arrays sized for N=32
     * (the largest test).  The same buffers are reused for all tests;
     * generate_twiddle_factors() overwrites them each call with the
     * correct N/2 entries for that run.
     * ---------------------------------------------------------------- */
    int max_n = 32;
    float *tr = (float *)malloc((max_n / 2) * sizeof(float));
    float *ti = (float *)malloc((max_n / 2) * sizeof(float));
    if (!tr || !ti) { fprintf(stderr, "malloc failed\n"); return 1; }

    /* Input: DC impulse [1, 0, 0, ..., 0]  →  FFT output: all 1+0i */
    float in_r[32] = {1.0f, 0};   /* {1.0f} sets [0]=1, rest are 0-initialised */
    float in_i[32] = {0};

    unsigned long long cyc[6];
    cyc[0] = run_test("\n=== ITERATIVE FFT (N=8) ===",
                      in_r, in_i, tr, ti,  8, 1);
    cyc[1] = run_test("\n=== RECURSIVE FFT (N=8) ===",
                      in_r, in_i, tr, ti,  8, 0);
    cyc[2] = run_test("\n=== ITERATIVE FFT (N=16) ===",
                      in_r, in_i, tr, ti, 16, 1);
    cyc[3] = run_test("\n=== RECURSIVE FFT (N=16) ===",
                      in_r, in_i, tr, ti, 16, 0);
    cyc[4] = run_test("\n=== ITERATIVE FFT (N=32) ===",
                      in_r, in_i, tr, ti, 32, 1);
    cyc[5] = run_test("\n=== RECURSIVE FFT (N=32) ===",
                      in_r, in_i, tr, ti, 32, 0);

    free(tr);
    free(ti);

    /* 5. Performance table */
    printf("\n=== PERFORMANCE COMPARISON ===\n");
    printf("%-12s %-6s %-15s %-20s\n",
           "Impl","N","Clock Cycles","Stack Memory (bytes)");
    printf("%-12s %-6s %-15s %-20s\n",
           "------------","------",
           "---------------","--------------------");
    int sizes[3] = {8, 16, 32};
    for (int i = 0; i < 3; i++) {
        printf("%-12s %-6d %-15llu %-20d\n",
               "Iterative", sizes[i], cyc[i*2],   iter_stack(sizes[i]));
        printf("%-12s %-6d %-15llu %-20d\n",
               "Recursive", sizes[i], cyc[i*2+1], recu_stack(sizes[i]));
    }

    printf("\n--- Stack Memory Derivation ---\n");
    printf("Iterative: 176 bytes (fixed, no recursion)\n");
    printf("  butterfly_iterative frame : 80 bytes (ra + s0-s8, 10 regs x 8)\n");
    printf("  bit_reverse_array frame   : 48 bytes (ra + s0-s4, 6 regs x 8)\n");
    printf("  log2_int (leaf call)      : ~16 bytes (ra only)\n");
    printf("  padding / call overhead   : ~32 bytes\n");
    printf("  Total (fixed depth=3)     : 176 bytes\n");
    printf("\nRecursive: 96 + 64*log2(N) bytes\n");
    printf("  64 = one butterfly_recursive_helper frame (ra + s0-s6, 8 regs x 8)\n");
    printf("  96 = bit_reverse_array(48) + log2_int(16) + top-level overhead(32)\n");
    printf("  N=8  -> %d bytes\n", recu_stack(8));
    printf("  N=16 -> %d bytes\n", recu_stack(16));
    printf("  N=32 -> %d bytes\n", recu_stack(32));

    return 0;
}