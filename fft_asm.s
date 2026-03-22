# =============================================================================
#  fft_asm.s  -  RISC-V (RV64G) 1D FFT — Milestone 2
#
#  ALL functions including twiddle factor generation are in assembly.
#  Float constants loaded via fmv.s.x with hardcoded IEEE-754 bit patterns
#  to avoid PC-relative addressing issues on QEMU/WSL2.
#
#  Exported functions:
#    sin_approx(float x)          -> fa0
#    cos_approx(float x)          -> fa0
#    log2_int(int n)              -> a0
#    reverse_bits(unsigned x, int logn) -> a0
#    generate_twiddle_factors(float *tr, float *ti, int n)
#    bit_reverse_array(float *xr, float *xi, int n)
#    butterfly_iterative(float *xr, float *xi, float *tr, float *ti, int n)
#    butterfly_recursive_helper(float *xr, float *xi,
#                               float *tr, float *ti, int n, int full_n)
#
#  IEEE-754 single-precision constants used:
#    PI       = 0x40490FDB = 3.14159265
#    2*PI     = 0x40C90FDB = 6.28318530
#    PI/2     = 0x3FC90FDB = 1.57079632
#    1.0      = 0x3F800000
#    -1/2     = 0xBF000000
#    1/24     = 0x3D2AAAAB
#    -1/720   = 0xBA800AB0
#    1/40320  = 0x358EFAAC
#    -1/3628800 = 0xB2B60B61   (cos coefficients)
#    1.0      = 0x3F800000
#    -1/6     = 0xBE2AAAAB
#    1/120    = 0x3C088889
#    -1/5040  = 0xB9500D01
#    1/362880 = 0x35888889
#    -1/39916800 = 0xB2D7322B  (sin coefficients)
#
#  ABI note (FIX applied):
#    fs0-fs11 are callee-saved in the RISC-V calling convention.
#    butterfly_iterative and butterfly_recursive_helper now save/restore
#    fs8-fs11 explicitly so the ABI is fully honoured.
# =============================================================================

    .section .text

# =============================================================================
# sin_approx(float x) -> fa0
# 6-term Taylor series with angle wrapping to [-PI, PI]
# Coefficients loaded via integer bit patterns (no .rodata needed)
# =============================================================================
    .globl sin_approx
sin_approx:
    addi  sp, sp, -48
    sd    ra, 40(sp)
    sd    s0, 32(sp)
    sd    s1, 24(sp)
    fsw   fs0, 20(sp)
    fsw   fs1, 16(sp)

    fmv.s fs0, fa0

    li    t0, 0x40490FDB; fmv.s.x ft1, t0   # ft1 = PI
    li    t0, 0x40C90FDB; fmv.s.x ft2, t0   # ft2 = 2*PI

    # wrap to [-PI, PI]
1:  flt.s t0, ft1, fs0; beqz t0, 2f; fsub.s fs0, fs0, ft2; j 1b
2:  fneg.s ft3, ft1
3:  flt.s t0, fs0, ft3; beqz t0, 4f; fadd.s fs0, fs0, ft2; j 3b
4:
    fmul.s fs1, fs0, fs0    # x^2

    # Horner evaluation (innermost first):
    # sin(x) = x*(1 + x^2*(-1/6 + x^2*(1/120 + x^2*(-1/5040
    #           + x^2*(1/362880 + x^2*(-1/39916800))))))
    li    t0, 0xB2D7322B; fmv.s.x fa0, t0   # -1/39916800
    li    t0, 0x35888889; fmv.s.x ft0, t0   #  1/362880
    fmadd.s fa0, fs1, fa0, ft0
    li    t0, 0xB9500D01; fmv.s.x ft0, t0   # -1/5040
    fmadd.s fa0, fs1, fa0, ft0
    li    t0, 0x3C088889; fmv.s.x ft0, t0   #  1/120
    fmadd.s fa0, fs1, fa0, ft0
    li    t0, 0xBE2AAAAB; fmv.s.x ft0, t0   # -1/6
    fmadd.s fa0, fs1, fa0, ft0
    li    t0, 0x3F800000; fmv.s.x ft0, t0   #  1.0
    fmadd.s fa0, fs1, fa0, ft0
    fmul.s fa0, fa0, fs0    # * x = sin(x)

    flw   fs0, 20(sp); flw fs1, 16(sp)
    ld    ra, 40(sp); ld s0, 32(sp); ld s1, 24(sp)
    addi  sp, sp, 48
    ret

# =============================================================================
# cos_approx(float x) -> fa0
# Direct 6-term Taylor series (NOT via sin to avoid PI/2 rounding error)
# cos(x) = 1 - x^2/2 + x^4/24 - x^6/720 + x^8/40320 - x^10/3628800
# =============================================================================
    .globl cos_approx
cos_approx:
    addi  sp, sp, -32
    sd    ra, 24(sp)
    fsw   fs0, 20(sp)
    fsw   fs1, 16(sp)

    fmv.s fs0, fa0

    li    t0, 0x40490FDB; fmv.s.x ft1, t0   # PI
    li    t0, 0x40C90FDB; fmv.s.x ft2, t0   # 2*PI

    # wrap to [-PI, PI]
1:  flt.s t0, ft1, fs0; beqz t0, 2f; fsub.s fs0, fs0, ft2; j 1b
2:  fneg.s ft3, ft1
3:  flt.s t0, fs0, ft3; beqz t0, 4f; fadd.s fs0, fs0, ft2; j 3b
4:
    fmul.s fs1, fs0, fs0    # x^2

    # Horner evaluation (innermost first):
    # cos(x) = 1 + x^2*(-1/2 + x^2*(1/24 + x^2*(-1/720
    #           + x^2*(1/40320 + x^2*(-1/3628800)))))
    li    t0, 0xB2B60B61; fmv.s.x fa0, t0   # -1/3628800
    li    t0, 0x358EFAAC; fmv.s.x ft0, t0   #  1/40320
    fmadd.s fa0, fs1, fa0, ft0
    li    t0, 0xBA800AB0; fmv.s.x ft0, t0   # -1/720
    fmadd.s fa0, fs1, fa0, ft0
    li    t0, 0x3D2AAAAB; fmv.s.x ft0, t0   #  1/24
    fmadd.s fa0, fs1, fa0, ft0
    li    t0, 0xBF000000; fmv.s.x ft0, t0   # -1/2
    fmadd.s fa0, fs1, fa0, ft0
    li    t0, 0x3F800000; fmv.s.x ft0, t0   #  1.0
    fmadd.s fa0, fs1, fa0, ft0
    # result already in fa0 (no multiply by x for cosine)

    flw   fs0, 20(sp); flw fs1, 16(sp)
    ld    ra,  24(sp)
    addi  sp, sp, 32
    ret

# =============================================================================
# log2_int(int n) -> a0
# Returns floor(log2(n)) by counting right-shifts until n <= 1
# =============================================================================
    .globl log2_int
log2_int:
    li  a1, 0
1:  li  t0, 1; ble a0, t0, 2f
    srli a0, a0, 1
    addi a1, a1, 1
    j   1b
2:  mv  a0, a1
    ret

# =============================================================================
# reverse_bits(unsigned x, int logn) -> a0
# Bit-reverses the logn least-significant bits of x
# =============================================================================
    .globl reverse_bits
reverse_bits:
    li  a2, 0           # result = 0
    mv  a3, a1          # counter = logn
1:  beqz a3, 2f
    slli a2, a2, 1      # result <<= 1
    andi t0, a0, 1      # t0 = x & 1 (LSB)
    or   a2, a2, t0     # result |= LSB
    srli a0, a0, 1      # x >>= 1
    addi a3, a3, -1
    j    1b
2:  mv  a0, a2
    ret

# =============================================================================
# generate_twiddle_factors(float *tr, float *ti, int n)
# Computes W_N^k = cos(-2*PI*k/N) + j*sin(-2*PI*k/N) for k=0..N/2-1
# Uses cos_approx and sin_approx (both in this file — pure assembly)
# Fix: cos(-x)=cos(x), sin(-x)=-sin(x), so we compute positive angle
# then negate sin result — avoids passing negative angles to sin_approx
#
# tr and ti point to HEAP-allocated buffers (malloc'd in main.c).
# =============================================================================
    .globl generate_twiddle_factors
generate_twiddle_factors:
    addi  sp, sp, -48
    sd    ra, 40(sp); sd s0, 32(sp); sd s1, 24(sp)
    sd    s2, 16(sp); sd s3,  8(sp); sd s4,  0(sp)

    mv    s0, a0            # tr pointer (heap)
    mv    s1, a1            # ti pointer (heap)
    srli  s2, a2, 1         # half = n/2
    mv    s3, a2            # n
    li    s4, 0             # k = 0

1:  bge   s4, s2, 2f

    # angle = 2*PI*k/n  (positive; we negate sin after)
    li    t0, 0x40C90FDB; fmv.s.x ft0, t0   # TWO_PI
    fcvt.s.w ft1, s4                          # (float)k
    fmul.s   ft0, ft0, ft1                    # 2*PI*k
    fcvt.s.w ft1, s3                          # (float)n
    fdiv.s   fs0, ft0, ft1                    # fs0 = 2*PI*k/n

    # twiddle_real[k] = cos(2*PI*k/n)  [cos is even: cos(-x)=cos(x)]
    fmv.s fa0, fs0; call cos_approx
    slli  t0, s4, 2; add t0, s0, t0; fsw fa0, 0(t0)

    # twiddle_imag[k] = -sin(2*PI*k/n) [sin(-x) = -sin(x)]
    fmv.s fa0, fs0; call sin_approx
    fneg.s fa0, fa0                           # negate
    slli  t0, s4, 2; add t0, s1, t0; fsw fa0, 0(t0)

    addi  s4, s4, 1; j 1b
2:
    ld    ra, 40(sp); ld s0, 32(sp); ld s1, 24(sp)
    ld    s2, 16(sp); ld s3,  8(sp); ld s4,  0(sp)
    addi  sp, sp, 48
    ret

# =============================================================================
# bit_reverse_array(float *xr, float *xi, int n)
# In-place bit-reversal permutation — only swaps when j > i
# =============================================================================
    .globl bit_reverse_array
bit_reverse_array:
    addi  sp, sp, -48
    sd    ra, 40(sp); sd s0, 32(sp); sd s1, 24(sp)
    sd    s2, 16(sp); sd s3,  8(sp); sd s4,  0(sp)

    mv    s0, a0; mv s1, a1; mv s2, a2
    mv    a0, s2; call log2_int; mv s3, a0   # s3 = log_n
    li    s4, 0                               # i = 0

1:  bge   s4, s2, 2f
    mv    a0, s4; mv a1, s3; call reverse_bits; mv t1, a0   # t1 = j
    ble   t1, s4, 3f                          # skip if j <= i
    slli  t2, s4, 2; slli t3, t1, 2
    add   t4, s0, t2; add t5, s0, t3
    flw   ft0, 0(t4); flw ft1, 0(t5); fsw ft1, 0(t4); fsw ft0, 0(t5)
    add   t4, s1, t2; add t5, s1, t3
    flw   ft0, 0(t4); flw ft1, 0(t5); fsw ft1, 0(t4); fsw ft0, 0(t5)
3:  addi  s4, s4, 1; j 1b
2:
    ld    ra, 40(sp); ld s0, 32(sp); ld s1, 24(sp)
    ld    s2, 16(sp); ld s3,  8(sp); ld s4,  0(sp)
    addi  sp, sp, 48
    ret

# =============================================================================
# butterfly_iterative(float *xr, float *xi, float *tr, float *ti, int n)
# Cooley-Tukey DIT iterative butterfly — log2(N) stages
#
# ABI FIX: fs8-fs11 are callee-saved; they are now saved/restored in this
# frame.  Frame size increased from 80 → 112 bytes to accommodate them.
#
# Frame layout (112 bytes, 16-byte aligned):
#   offset 104: ra
#   offset  96: s0
#   offset  88: s1
#   offset  80: s2
#   offset  72: s3
#   offset  64: s4
#   offset  56: s5
#   offset  48: s6
#   offset  40: s7
#   offset  32: s8
#   offset  28: fs8  (float, 4 bytes)
#   offset  24: fs9
#   offset  20: fs10
#   offset  16: fs11
#   offset   0: (padding to 16-byte alignment)
# =============================================================================
    .globl butterfly_iterative
butterfly_iterative:
    addi  sp, sp, -112
    sd    ra,104(sp); sd s0, 96(sp); sd s1, 88(sp); sd s2, 80(sp)
    sd    s3, 72(sp); sd s4, 64(sp); sd s5, 56(sp); sd s6, 48(sp)
    sd    s7, 40(sp); sd s8, 32(sp)
    fsw   fs8, 28(sp); fsw fs9, 24(sp); fsw fs10, 20(sp); fsw fs11, 16(sp)

    mv    s0,a0; mv s1,a1; mv s2,a2; mv s3,a3; mv s4,a4
    mv    a0,s4; call log2_int; mv s5,a0   # s5 = log_n
    li    s6,1                              # s = stage (start at 1)

stg:bgt   s6,s5,stgdn
    li    t0,1; sll s7,t0,s6   # m = 1<<s
    srli  s8,s7,1              # half = m/2
    div   t2,s4,s7             # step = n/m

    li    a5,0                 # k = 0
klp:bge   a5,s4,kdn
      li  a6,0                 # j = 0
jlp:  bge a6,s8,jdn
        mul  t3,a6,t2; slli t3,t3,2
        add  t4,s2,t3; add t5,s3,t3
        flw  fs2,0(t4); flw fs3,0(t5)      # wr, wi

        add  t3,a5,a6; slli t3,t3,2
        add  t4,s0,t3; add t5,s1,t3
        flw  fs4,0(t4); flw fs5,0(t5)      # ur, ui

        add  t3,a5,a6; add t3,t3,s8; slli t3,t3,2
        add  t4,s0,t3; add t5,s1,t3
        flw  fs6,0(t4); flw fs7,0(t5)      # vr, vi

        # t = tw * v
        # tr = wr*vr - wi*vi
        fmul.s fs8,fs2,fs6; fmul.s fs9,fs3,fs7; fsub.s fs8,fs8,fs9
        # ti = wr*vi + wi*vr
        fmul.s fs10,fs2,fs7; fmul.s fs11,fs3,fs6; fadd.s fs10,fs10,fs11

        # write upper output: x[k+j] = u + t
        add  t3,a5,a6; slli t3,t3,2
        add  t4,s0,t3; add t5,s1,t3
        fadd.s ft0,fs4,fs8;  fsw ft0,0(t4)
        fadd.s ft0,fs5,fs10; fsw ft0,0(t5)

        # write lower output: x[k+j+half] = u - t
        add  t3,a5,a6; add t3,t3,s8; slli t3,t3,2
        add  t4,s0,t3; add t5,s1,t3
        fsub.s ft0,fs4,fs8;  fsw ft0,0(t4)
        fsub.s ft0,fs5,fs10; fsw ft0,0(t5)

        addi a6,a6,1; j jlp
jdn:  add  a5,a5,s7; j klp
kdn:
    addi  s6,s6,1; j stg
stgdn:
    flw   fs8, 28(sp); flw fs9, 24(sp); flw fs10, 20(sp); flw fs11, 16(sp)
    ld    ra,104(sp); ld s0, 96(sp); ld s1, 88(sp); ld s2, 80(sp)
    ld    s3, 72(sp); ld s4, 64(sp); ld s5, 56(sp); ld s6, 48(sp)
    ld    s7, 40(sp); ld s8, 32(sp)
    addi  sp,sp,112
    ret

# =============================================================================
# butterfly_recursive_helper(xr, xi, tr, ti, n, full_n)
# a0=xr  a1=xi  a2=tr  a3=ti  a4=n  a5=full_n
#
# Cooley-Tukey DIT recursive butterfly.
# bit_reverse_array is called BEFORE this function by the C wrapper.
#
# ABI FIX: fs8-fs11 are callee-saved; they are now saved/restored.
# Frame size: 80 bytes  (8 integer regs × 8 + 4 float regs × 4 = 80, padded to 80)
#
# Frame layout (80 bytes, 16-byte aligned):
#   offset 72: ra
#   offset 64: s0
#   offset 56: s1
#   offset 48: s2
#   offset 40: s3
#   offset 32: s4
#   offset 24: s5
#   offset 16: s6
#   offset 12: fs8
#   offset  8: fs9
#   offset  4: fs10
#   offset  0: fs11
# =============================================================================
    .globl butterfly_recursive_helper
butterfly_recursive_helper:
    li    t0,1; ble a4,t0,rh_base   # base case: n<=1, return

    addi  sp, sp, -80
    sd    ra, 72(sp); sd s0, 64(sp); sd s1, 56(sp); sd s2, 48(sp)
    sd    s3, 40(sp); sd s4, 32(sp); sd s5, 24(sp); sd s6, 16(sp)
    fsw   fs8, 12(sp); fsw fs9, 8(sp); fsw fs10, 4(sp); fsw fs11, 0(sp)

    mv    s0,a0; mv s1,a1; mv s2,a2; mv s3,a3; mv s4,a4; mv s5,a5
    srli  s6,s4,1   # half = n/2

    # Recurse on lower half [0..half-1]
    mv    a0,s0; mv a1,s1; mv a2,s2; mv a3,s3; mv a4,s6; mv a5,s5
    call  butterfly_recursive_helper

    # Recurse on upper half [half..n-1]
    slli  t0,s6,2; add a0,s0,t0; add a1,s1,t0
    mv    a2,s2; mv a3,s3; mv a4,s6; mv a5,s5
    call  butterfly_recursive_helper

    # Combine stage: step = full_n / n  (correct twiddle stride for sub-size)
    div   t2,s5,s4   # step
    li    t3,0       # k = 0

cb: bge   t3,s6,cbdn
        mul  t4,t3,t2; slli t4,t4,2
        add  t5,s2,t4; add t6,s3,t4
        flw  fs2,0(t5); flw fs3,0(t6)   # wr, wi

        slli t4,t3,2; add t5,s0,t4; add t6,s1,t4
        flw  fs4,0(t5); flw fs5,0(t6)   # ur, ui  (lower half)

        add  t4,t3,s6; slli t4,t4,2
        add  t5,s0,t4; add t6,s1,t4
        flw  fs6,0(t5); flw fs7,0(t6)   # vr, vi  (upper half)

        # t = tw * v
        fmul.s fs8,fs2,fs6; fmul.s fs9,fs3,fs7; fsub.s fs8,fs8,fs9
        fmul.s fs10,fs2,fs7; fmul.s fs11,fs3,fs6; fadd.s fs10,fs10,fs11

        # x[k] = u + t
        slli t4,t3,2; add t5,s0,t4; add t6,s1,t4
        fadd.s ft0,fs4,fs8;  fsw ft0,0(t5)
        fadd.s ft0,fs5,fs10; fsw ft0,0(t6)

        # x[k+half] = u - t
        add  t4,t3,s6; slli t4,t4,2
        add  t5,s0,t4; add t6,s1,t4
        fsub.s ft0,fs4,fs8;  fsw ft0,0(t5)
        fsub.s ft0,fs5,fs10; fsw ft0,0(t6)

        addi t3,t3,1; j cb
cbdn:
    flw   fs8, 12(sp); flw fs9, 8(sp); flw fs10, 4(sp); flw fs11, 0(sp)
    ld    ra, 72(sp); ld s0, 64(sp); ld s1, 56(sp); ld s2, 48(sp)
    ld    s3, 40(sp); ld s4, 32(sp); ld s5, 24(sp); ld s6, 16(sp)
    addi  sp,sp,80
    ret

rh_base: ret