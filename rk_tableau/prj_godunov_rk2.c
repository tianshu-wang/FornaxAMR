#include "prj.h"

/*
 * prj_godunov_rk2 -- an additive (IMEX) Butcher tableau that reproduces the
 * hardwired TIME_INTEGRATION=RK2 path bit-for-bit (up to minor differences
 * that the RK2 code path glosses over, e.g. the lapse not being re-evaluated
 * inside the implicit radiation solve).
 *
 * TIME_INTEGRATION=RK2 is an operator split: in each of the two RK2 stages it
 * first does an explicit transport advance (flux divergence + non-stiff
 * sources, "F") and then applies the radiation source as a backward-Euler
 * implicit solve ("G", prj_rad_energy_update & co., u_out = u_in + dt G(u_out)).
 * Stage 1 uses the full dt for the implicit solve, stage 2 uses 0.5*dt.
 *
 * Writing du/dt = F(u) + G(u), that split is exactly the 3-stage additive RK
 * below, with a_ex/b_ex the explicit (transport) coefficients and
 * a_im/b_im the implicit (radiation) coefficients:
 *
 *   stage 1:  U1 = u^n                         (evaluate F(u^n); no implicit solve)
 *   stage 2:  U2 = u^n + dt F(u^n) + dt G(U2)  (= RK2 stage 1, implicit dt)
 *   stage 3:  U3 = u^n + dt(0.5 F(U1)+0.5 F(U2)) + dt(0.5 G(U2)+0.5 G(U3))
 *                                              (= RK2 stage 2, implicit 0.5 dt)
 *
 * The weights equal the last stage rows (b_ex = a_ex[2][:], b_im = a_im[2][:]),
 * so the method is stiffly accurate / FSAL: u^{n+1} = U3, matching the RK2
 * code whose step output is its stage-2 result.
 *
 * Order note: the explicit part is 2nd-order (Heun) and each implicit solve is
 * backward Euler, but the additive coupling conditions are not met
 * (sum b_im_i c_ex_i = 1 != 1/2), so -- like any Lie/Godunov operator split --
 * the combined scheme is only 1st-order accurate. It is provided as a faithful
 * tableau encoding of the existing RK2 behavior, not as a 2nd-order IMEX method.
 *
 * Explicit A:            Implicit A:
 *   [ 0    0    0 ]        [ 0    0    0   ]
 *   [ 1    0    0 ]        [ 0    1    0   ]
 *   [ 0.5  0.5  0 ]        [ 0    0.5  0.5 ]
 * b_ex = [0.5, 0.5, 0]   b_im = [0, 0.5, 0.5]
 */

const double prj_godunov_rk2_a_ex[3 * 3] = {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    0.5, 0.5, 0.0
};

const double prj_godunov_rk2_a_im[3 * 3] = {
    0.0, 0.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.5, 0.5
};

const double prj_godunov_rk2_b_ex[3] = {
    0.5,
    0.5,
    0.0
};

const double prj_godunov_rk2_b_im[3] = {
    0.0,
    0.5,
    0.5
};

const prj_timeint_imex_tableau prj_godunov_rk2 = {
    3,
    prj_godunov_rk2_a_ex,
    prj_godunov_rk2_a_im,
    prj_godunov_rk2_b_ex,
    prj_godunov_rk2_b_im,
    1.0
};
