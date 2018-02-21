/*
 * helikopter_4.c
 *
 * Code generation for model "helikopter_4".
 *
 * Model version              : 1.178
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Tue Feb 20 16:45:31 2018
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "helikopter_4.h"
#include "helikopter_4_private.h"
#include "helikopter_4_dt.h"

/* Block signals (auto storage) */
B_helikopter_4_T helikopter_4_B;

/* Continuous states */
X_helikopter_4_T helikopter_4_X;

/* Block states (auto storage) */
DW_helikopter_4_T helikopter_4_DW;

/* Real-time model */
RT_MODEL_helikopter_4_T helikopter_4_M_;
RT_MODEL_helikopter_4_T *const helikopter_4_M = &helikopter_4_M_;

/*
 * This function updates continuous states using the ODE1 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE1_IntgData *id = (ODE1_IntgData *)rtsiGetSolverData(si);
  real_T *f0 = id->f[0];
  int_T i;
  int_T nXc = 4;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  helikopter_4_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helikopter_4_output(void)
{
  /* local block i/o variables */
  real_T rtb_FromWorkspace[2];
  real_T rtb_Backgain;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T *lastU;
  real_T rtb_Gain1_idx_4;
  real_T rtb_Gain1_idx_5;
  if (rtmIsMajorTimeStep(helikopter_4_M)) {
    /* set solver stop time */
    if (!(helikopter_4_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helikopter_4_M->solverInfo,
                            ((helikopter_4_M->Timing.clockTickH0 + 1) *
        helikopter_4_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helikopter_4_M->solverInfo,
                            ((helikopter_4_M->Timing.clockTick0 + 1) *
        helikopter_4_M->Timing.stepSize0 + helikopter_4_M->Timing.clockTickH0 *
        helikopter_4_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helikopter_4_M)) {
    helikopter_4_M->Timing.t[0] = rtsiGetT(&helikopter_4_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helikopter_4_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S5>/HIL Read Encoder Timebase' */

    /* S-Function Block: helikopter_4/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder(helikopter_4_DW.HILReadEncoderTimebase_Task,
        1, &helikopter_4_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helikopter_4_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helikopter_4_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          helikopter_4_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *) helikopter_4_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) helikopter_4_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = helikopter_4_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = helikopter_4_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[80]) {
      currTimeIndex = 79;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helikopter_4_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 2; ++elIdx) {
              (&rtb_FromWorkspace[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 81;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 2; ++elIdx) {
              (&rtb_FromWorkspace[0])[elIdx] = pDataValues[currTimeIndex + 1];
              pDataValues += 81;
            }
          }
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;

        {
          int_T elIdx;
          for (elIdx = 0; elIdx < 2; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_FromWorkspace[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1,
              f2);
            pDataValues += 81;
          }
        }
      }
    }
  }

  /* Sum: '<Root>/Sum4' */
  helikopter_4_B.u[0] = rtb_FromWorkspace[0] - 0.0;
  helikopter_4_B.u[1] = rtb_FromWorkspace[1] - 0.0;
  if (rtmIsMajorTimeStep(helikopter_4_M)) {
    /* Gain: '<S5>/Travel: Count to rad' */
    helikopter_4_B.TravelCounttorad = helikopter_4_P.TravelCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o1;

    /* Gain: '<S12>/Gain' */
    helikopter_4_B.Gain = helikopter_4_P.Gain_Gain *
      helikopter_4_B.TravelCounttorad;

    /* Sum: '<Root>/Sum5' incorporates:
     *  Constant: '<Root>/Constant'
     */
    helikopter_4_B.travel = helikopter_4_P.Constant_Value + helikopter_4_B.Gain;
  }

  /* TransferFcn: '<S5>/Travel: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helikopter_4_P.TravelTransferFcn_C *
    helikopter_4_X.TravelTransferFcn_CSTATE;
  rtb_Backgain += helikopter_4_P.TravelTransferFcn_D *
    helikopter_4_B.TravelCounttorad;

  /* Gain: '<S13>/Gain' */
  helikopter_4_B.Gain_d = helikopter_4_P.Gain_Gain_l * rtb_Backgain;
  if (rtmIsMajorTimeStep(helikopter_4_M)) {
    /* Gain: '<S5>/Pitch: Count to rad' */
    helikopter_4_B.PitchCounttorad = helikopter_4_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S9>/Gain' */
    helikopter_4_B.Gain_i = helikopter_4_P.Gain_Gain_a *
      helikopter_4_B.PitchCounttorad;
  }

  /* TransferFcn: '<S5>/Pitch: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helikopter_4_P.PitchTransferFcn_C *
    helikopter_4_X.PitchTransferFcn_CSTATE;
  rtb_Backgain += helikopter_4_P.PitchTransferFcn_D *
    helikopter_4_B.PitchCounttorad;

  /* Gain: '<S10>/Gain' */
  helikopter_4_B.Gain_b = helikopter_4_P.Gain_Gain_ae * rtb_Backgain;
  if (rtmIsMajorTimeStep(helikopter_4_M)) {
    /* Gain: '<S5>/Elevation: Count to rad' */
    helikopter_4_B.ElevationCounttorad = helikopter_4_P.ElevationCounttorad_Gain
      * rtb_HILReadEncoderTimebase_o3;

    /* Gain: '<S7>/Gain' */
    helikopter_4_B.Gain_e = helikopter_4_P.Gain_Gain_lv *
      helikopter_4_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helikopter_4_B.Sum = helikopter_4_B.Gain_e +
      helikopter_4_P.elavation_offsetdeg_Value;
  }

  /* TransferFcn: '<S5>/Elevation: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helikopter_4_P.ElevationTransferFcn_C *
    helikopter_4_X.ElevationTransferFcn_CSTATE;
  rtb_Backgain += helikopter_4_P.ElevationTransferFcn_D *
    helikopter_4_B.ElevationCounttorad;

  /* Gain: '<S8>/Gain' */
  helikopter_4_B.Gain_dg = helikopter_4_P.Gain_Gain_n * rtb_Backgain;
  if (rtmIsMajorTimeStep(helikopter_4_M)) {
    /* SignalConversion: '<Root>/TmpSignal ConversionAtTo WorkspaceInport1' */
    helikopter_4_B.TmpSignalConversionAtToWorkspac[0] = helikopter_4_B.u[0];
    helikopter_4_B.TmpSignalConversionAtToWorkspac[1] = helikopter_4_B.u[1];
    helikopter_4_B.TmpSignalConversionAtToWorkspac[2] = helikopter_4_B.travel;
    helikopter_4_B.TmpSignalConversionAtToWorkspac[3] = helikopter_4_B.Gain_d;
    helikopter_4_B.TmpSignalConversionAtToWorkspac[4] = helikopter_4_B.Gain_i;
    helikopter_4_B.TmpSignalConversionAtToWorkspac[5] = helikopter_4_B.Gain_b;
    helikopter_4_B.TmpSignalConversionAtToWorkspac[6] = helikopter_4_B.Sum;
    helikopter_4_B.TmpSignalConversionAtToWorkspac[7] = helikopter_4_B.Gain_dg;
  }

  /* Gain: '<S2>/Gain1' */
  rtb_Gain1_idx_4 = helikopter_4_P.Gain1_Gain * helikopter_4_B.Sum;
  rtb_Gain1_idx_5 = helikopter_4_P.Gain1_Gain * helikopter_4_B.Gain_dg;

  /* Sum: '<Root>/Sum1' incorporates:
   *  Constant: '<Root>/Vd_bias'
   *  Gain: '<S2>/Gain1'
   *  Gain: '<S6>/K_pd'
   *  Gain: '<S6>/K_pp'
   *  Sum: '<S6>/Sum2'
   *  Sum: '<S6>/Sum3'
   */
  helikopter_4_B.Sum1 = ((helikopter_4_B.u[0] - helikopter_4_P.Gain1_Gain *
    helikopter_4_B.Gain_i) * helikopter_4_P.K_pp - helikopter_4_P.Gain1_Gain *
    helikopter_4_B.Gain_b * helikopter_4_P.K_pd) + helikopter_4_P.Vd_ff;
  if (rtmIsMajorTimeStep(helikopter_4_M)) {
  }

  /* Integrator: '<S4>/Integrator'
   *
   * Regarding '<S4>/Integrator':
   *  Limited Integrator
   */
  if (helikopter_4_X.Integrator_CSTATE >= helikopter_4_P.Integrator_UpperSat ) {
    helikopter_4_X.Integrator_CSTATE = helikopter_4_P.Integrator_UpperSat;
  } else if (helikopter_4_X.Integrator_CSTATE <=
             (helikopter_4_P.Integrator_LowerSat) ) {
    helikopter_4_X.Integrator_CSTATE = (helikopter_4_P.Integrator_LowerSat);
  }

  rtb_Backgain = helikopter_4_X.Integrator_CSTATE;

  /* Sum: '<S4>/Sum' */
  rtb_Gain1_idx_4 = helikopter_4_B.u[1] - rtb_Gain1_idx_4;

  /* Sum: '<Root>/Sum2' incorporates:
   *  Constant: '<Root>/Vs_bias'
   *  Gain: '<S4>/K_ed'
   *  Gain: '<S4>/K_ep'
   *  Sum: '<S4>/Sum1'
   */
  helikopter_4_B.Sum2 = ((helikopter_4_P.K_ep * rtb_Gain1_idx_4 + rtb_Backgain)
    - helikopter_4_P.K_ed * rtb_Gain1_idx_5) + helikopter_4_P.Vs_ff;
  if (rtmIsMajorTimeStep(helikopter_4_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Backgain = (helikopter_4_B.Sum2 - helikopter_4_B.Sum1) *
    helikopter_4_P.Backgain_Gain;

  /* Gain: '<S4>/K_ei' */
  helikopter_4_B.K_ei = helikopter_4_P.K_ei * rtb_Gain1_idx_4;
  if (rtmIsMajorTimeStep(helikopter_4_M)) {
  }

  /* Derivative: '<S5>/Derivative' */
  if ((helikopter_4_DW.TimeStampA >= helikopter_4_M->Timing.t[0]) &&
      (helikopter_4_DW.TimeStampB >= helikopter_4_M->Timing.t[0])) {
    rtb_Gain1_idx_4 = 0.0;
  } else {
    rtb_Gain1_idx_4 = helikopter_4_DW.TimeStampA;
    lastU = &helikopter_4_DW.LastUAtTimeA;
    if (helikopter_4_DW.TimeStampA < helikopter_4_DW.TimeStampB) {
      if (helikopter_4_DW.TimeStampB < helikopter_4_M->Timing.t[0]) {
        rtb_Gain1_idx_4 = helikopter_4_DW.TimeStampB;
        lastU = &helikopter_4_DW.LastUAtTimeB;
      }
    } else {
      if (helikopter_4_DW.TimeStampA >= helikopter_4_M->Timing.t[0]) {
        rtb_Gain1_idx_4 = helikopter_4_DW.TimeStampB;
        lastU = &helikopter_4_DW.LastUAtTimeB;
      }
    }

    rtb_Gain1_idx_4 = (helikopter_4_B.PitchCounttorad - *lastU) /
      (helikopter_4_M->Timing.t[0] - rtb_Gain1_idx_4);
  }

  /* End of Derivative: '<S5>/Derivative' */

  /* Gain: '<S11>/Gain' */
  helikopter_4_B.Gain_l = helikopter_4_P.Gain_Gain_a1 * rtb_Gain1_idx_4;
  if (rtmIsMajorTimeStep(helikopter_4_M)) {
  }

  /* Saturate: '<S5>/Back motor: Saturation' */
  if (rtb_Backgain > helikopter_4_P.BackmotorSaturation_UpperSat) {
    helikopter_4_B.BackmotorSaturation =
      helikopter_4_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Backgain < helikopter_4_P.BackmotorSaturation_LowerSat) {
    helikopter_4_B.BackmotorSaturation =
      helikopter_4_P.BackmotorSaturation_LowerSat;
  } else {
    helikopter_4_B.BackmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S5>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helikopter_4_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Gain1_idx_4 = (helikopter_4_B.Sum1 + helikopter_4_B.Sum2) *
    helikopter_4_P.Frontgain_Gain;

  /* Saturate: '<S5>/Front motor: Saturation' */
  if (rtb_Gain1_idx_4 > helikopter_4_P.FrontmotorSaturation_UpperSat) {
    helikopter_4_B.FrontmotorSaturation =
      helikopter_4_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Gain1_idx_4 < helikopter_4_P.FrontmotorSaturation_LowerSat) {
    helikopter_4_B.FrontmotorSaturation =
      helikopter_4_P.FrontmotorSaturation_LowerSat;
  } else {
    helikopter_4_B.FrontmotorSaturation = rtb_Gain1_idx_4;
  }

  /* End of Saturate: '<S5>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helikopter_4_M)) {
    /* S-Function (hil_write_analog_block): '<S5>/HIL Write Analog' */

    /* S-Function Block: helikopter_4/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helikopter_4_DW.HILWriteAnalog_Buffer[0] =
        helikopter_4_B.FrontmotorSaturation;
      helikopter_4_DW.HILWriteAnalog_Buffer[1] =
        helikopter_4_B.BackmotorSaturation;
      result = hil_write_analog(helikopter_4_DW.HILInitialize_Card,
        helikopter_4_P.HILWriteAnalog_channels, 2,
        &helikopter_4_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void helikopter_4_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S5>/Derivative' */
  if (helikopter_4_DW.TimeStampA == (rtInf)) {
    helikopter_4_DW.TimeStampA = helikopter_4_M->Timing.t[0];
    lastU = &helikopter_4_DW.LastUAtTimeA;
  } else if (helikopter_4_DW.TimeStampB == (rtInf)) {
    helikopter_4_DW.TimeStampB = helikopter_4_M->Timing.t[0];
    lastU = &helikopter_4_DW.LastUAtTimeB;
  } else if (helikopter_4_DW.TimeStampA < helikopter_4_DW.TimeStampB) {
    helikopter_4_DW.TimeStampA = helikopter_4_M->Timing.t[0];
    lastU = &helikopter_4_DW.LastUAtTimeA;
  } else {
    helikopter_4_DW.TimeStampB = helikopter_4_M->Timing.t[0];
    lastU = &helikopter_4_DW.LastUAtTimeB;
  }

  *lastU = helikopter_4_B.PitchCounttorad;

  /* End of Update for Derivative: '<S5>/Derivative' */
  if (rtmIsMajorTimeStep(helikopter_4_M)) {
    rt_ertODEUpdateContinuousStates(&helikopter_4_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++helikopter_4_M->Timing.clockTick0)) {
    ++helikopter_4_M->Timing.clockTickH0;
  }

  helikopter_4_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helikopter_4_M->solverInfo);

  {
    /* Update absolute timer for sample time: [0.002s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++helikopter_4_M->Timing.clockTick1)) {
      ++helikopter_4_M->Timing.clockTickH1;
    }

    helikopter_4_M->Timing.t[1] = helikopter_4_M->Timing.clockTick1 *
      helikopter_4_M->Timing.stepSize1 + helikopter_4_M->Timing.clockTickH1 *
      helikopter_4_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void helikopter_4_derivatives(void)
{
  XDot_helikopter_4_T *_rtXdot;
  _rtXdot = ((XDot_helikopter_4_T *) helikopter_4_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S5>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += helikopter_4_P.TravelTransferFcn_A *
    helikopter_4_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helikopter_4_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S5>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += helikopter_4_P.PitchTransferFcn_A *
    helikopter_4_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helikopter_4_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S5>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE += helikopter_4_P.ElevationTransferFcn_A *
    helikopter_4_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += helikopter_4_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S4>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helikopter_4_X.Integrator_CSTATE <=
            (helikopter_4_P.Integrator_LowerSat) );
    usat = ( helikopter_4_X.Integrator_CSTATE >=
            helikopter_4_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helikopter_4_B.K_ei > 0)) ||
        (usat && (helikopter_4_B.K_ei < 0)) ) {
      ((XDot_helikopter_4_T *) helikopter_4_M->ModelData.derivs)
        ->Integrator_CSTATE = helikopter_4_B.K_ei;
    } else {
      /* in saturation */
      ((XDot_helikopter_4_T *) helikopter_4_M->ModelData.derivs)
        ->Integrator_CSTATE = 0.0;
    }
  }
}

/* Model initialize function */
void helikopter_4_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopter_4/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &helikopter_4_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options(helikopter_4_DW.HILInitialize_Card,
      "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helikopter_4_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
      return;
    }

    if ((helikopter_4_P.HILInitialize_set_analog_input_ && !is_switching) ||
        (helikopter_4_P.HILInitialize_set_analog_inpu_m && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helikopter_4_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = helikopter_4_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helikopter_4_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = helikopter_4_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges(helikopter_4_DW.HILInitialize_Card,
        helikopter_4_P.HILInitialize_analog_input_chan, 8U,
        &helikopter_4_DW.HILInitialize_AIMinimums[0],
        &helikopter_4_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter_4_P.HILInitialize_set_analog_output && !is_switching) ||
        (helikopter_4_P.HILInitialize_set_analog_outp_b && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &helikopter_4_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = helikopter_4_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &helikopter_4_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = helikopter_4_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges(helikopter_4_DW.HILInitialize_Card,
        helikopter_4_P.HILInitialize_analog_output_cha, 8U,
        &helikopter_4_DW.HILInitialize_AOMinimums[0],
        &helikopter_4_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter_4_P.HILInitialize_set_analog_outp_e && !is_switching) ||
        (helikopter_4_P.HILInitialize_set_analog_outp_j && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helikopter_4_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helikopter_4_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(helikopter_4_DW.HILInitialize_Card,
        helikopter_4_P.HILInitialize_analog_output_cha, 8U,
        &helikopter_4_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
        return;
      }
    }

    if (helikopter_4_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helikopter_4_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helikopter_4_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helikopter_4_DW.HILInitialize_Card,
         helikopter_4_P.HILInitialize_analog_output_cha, 8U,
         &helikopter_4_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter_4_P.HILInitialize_set_encoder_param && !is_switching) ||
        (helikopter_4_P.HILInitialize_set_encoder_par_m && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helikopter_4_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = helikopter_4_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (helikopter_4_DW.HILInitialize_Card,
         helikopter_4_P.HILInitialize_encoder_channels, 8U,
         (t_encoder_quadrature_mode *)
         &helikopter_4_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter_4_P.HILInitialize_set_encoder_count && !is_switching) ||
        (helikopter_4_P.HILInitialize_set_encoder_cou_k && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helikopter_4_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] =
            helikopter_4_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts(helikopter_4_DW.HILInitialize_Card,
        helikopter_4_P.HILInitialize_encoder_channels, 8U,
        &helikopter_4_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter_4_P.HILInitialize_set_pwm_params_at && !is_switching) ||
        (helikopter_4_P.HILInitialize_set_pwm_params__f && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues = &helikopter_4_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helikopter_4_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(helikopter_4_DW.HILInitialize_Card,
        helikopter_4_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &helikopter_4_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          helikopter_4_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues = &helikopter_4_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            helikopter_4_DW.HILInitialize_POSortedChans[num_duty_cycle_modes] =
              p_HILInitialize_pwm_channels[i1];
            helikopter_4_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes] =
              helikopter_4_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            helikopter_4_DW.HILInitialize_POSortedChans[7U - num_frequency_modes]
              = p_HILInitialize_pwm_channels[i1];
            helikopter_4_DW.HILInitialize_POSortedFreqs[7U - num_frequency_modes]
              = helikopter_4_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(helikopter_4_DW.HILInitialize_Card,
          &helikopter_4_DW.HILInitialize_POSortedChans[0], num_duty_cycle_modes,
          &helikopter_4_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(helikopter_4_DW.HILInitialize_Card,
          &helikopter_4_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &helikopter_4_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues = &helikopter_4_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helikopter_4_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &helikopter_4_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = helikopter_4_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helikopter_4_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = helikopter_4_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration(helikopter_4_DW.HILInitialize_Card,
        helikopter_4_P.HILInitialize_pwm_channels, 8U,
        (t_pwm_configuration *) &helikopter_4_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &helikopter_4_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &helikopter_4_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs = &helikopter_4_DW.HILInitialize_POSortedFreqs
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = helikopter_4_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helikopter_4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helikopter_4_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband(helikopter_4_DW.HILInitialize_Card,
        helikopter_4_P.HILInitialize_pwm_channels, 8U,
        &helikopter_4_DW.HILInitialize_POSortedFreqs[0],
        &helikopter_4_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter_4_P.HILInitialize_set_pwm_outputs_a && !is_switching) ||
        (helikopter_4_P.HILInitialize_set_pwm_outputs_g && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helikopter_4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helikopter_4_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(helikopter_4_DW.HILInitialize_Card,
        helikopter_4_P.HILInitialize_pwm_channels, 8U,
        &helikopter_4_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
        return;
      }
    }

    if (helikopter_4_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues = &helikopter_4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helikopter_4_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helikopter_4_DW.HILInitialize_Card,
         helikopter_4_P.HILInitialize_pwm_channels, 8U,
         &helikopter_4_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S5>/HIL Read Encoder Timebase' */

  /* S-Function Block: helikopter_4/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader(helikopter_4_DW.HILInitialize_Card,
      helikopter_4_P.HILReadEncoderTimebase_samples_,
      helikopter_4_P.HILReadEncoderTimebase_channels, 3,
      &helikopter_4_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
    }
  }

  /* Start for FromWorkspace: '<Root>/From Workspace' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -2.6632718234239436E-22, -2.0146511265741569E-21, -3.1467187463805931E-22,
      6.0841737138787793E-21, 4.49852324367517E-21, 7.7677349374864188E-22,
      -2.4057320217562152E-20, 4.7193514714667247E-20, -4.9609022928582351E-20,
      3.1867168382832372E-20, -1.2303773938971594E-20, 3.3858072023778312E-20,
      1.358351228834184E-20, 5.7606901729490319E-21, -1.4163313349070169E-21,
      1.8694578184639294E-25, 3.406105774720244E-22, 1.3210217356464342E-21,
      -6.23909369468398E-22, 1.2138333487860573E-22, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -2.1364728032152258E-5, -2.0486130563977741E-5, -1.9710905013766882E-5,
      -1.9039145512653811E-5, -1.847032735596222E-5, -1.8003592745635002E-5,
      -1.763797740202636E-5, -1.7372586967670028E-5, -1.7206729587718606E-5,
      -1.7140010670996749E-5, -1.7172393686528524E-5, -1.7304229342321646E-5,
      -1.7536253400417833E-5, -1.7869554956810284E-5, -1.8305511234713053E-5,
      -1.8845688070816049E-5, -1.9491699625800659E-5, -2.0245024293212502E-5,
      -2.1106765347656437E-5, -2.2077350556028696E-5, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helikopter_4_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    helikopter_4_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    helikopter_4_DW.FromWorkspace_IWORK.PrevIndex = 0;
  }

  /* InitializeConditions for TransferFcn: '<S5>/Travel: Transfer Fcn' */
  helikopter_4_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S5>/Pitch: Transfer Fcn' */
  helikopter_4_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S5>/Elevation: Transfer Fcn' */
  helikopter_4_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S4>/Integrator' */
  helikopter_4_X.Integrator_CSTATE = helikopter_4_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S5>/Derivative' */
  helikopter_4_DW.TimeStampA = (rtInf);
  helikopter_4_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void helikopter_4_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopter_4/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helikopter_4_DW.HILInitialize_Card);
    hil_monitor_stop_all(helikopter_4_DW.HILInitialize_Card);
    is_switching = false;
    if ((helikopter_4_P.HILInitialize_set_analog_out_ex && !is_switching) ||
        (helikopter_4_P.HILInitialize_set_analog_outp_c && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helikopter_4_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helikopter_4_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((helikopter_4_P.HILInitialize_set_pwm_output_ap && !is_switching) ||
        (helikopter_4_P.HILInitialize_set_pwm_outputs_p && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helikopter_4_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helikopter_4_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helikopter_4_DW.HILInitialize_Card
                         , helikopter_4_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , helikopter_4_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helikopter_4_DW.HILInitialize_AOVoltages[0]
                         , &helikopter_4_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(helikopter_4_DW.HILInitialize_Card,
            helikopter_4_P.HILInitialize_analog_output_cha,
            num_final_analog_outputs, &helikopter_4_DW.HILInitialize_AOVoltages
            [0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(helikopter_4_DW.HILInitialize_Card,
            helikopter_4_P.HILInitialize_pwm_channels, num_final_pwm_outputs,
            &helikopter_4_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helikopter_4_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helikopter_4_DW.HILInitialize_Card);
    hil_monitor_delete_all(helikopter_4_DW.HILInitialize_Card);
    hil_close(helikopter_4_DW.HILInitialize_Card);
    helikopter_4_DW.HILInitialize_Card = NULL;
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  helikopter_4_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  helikopter_4_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  helikopter_4_initialize();
}

void MdlTerminate(void)
{
  helikopter_4_terminate();
}

/* Registration function */
RT_MODEL_helikopter_4_T *helikopter_4(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helikopter_4_P.Integrator_UpperSat = rtInf;
  helikopter_4_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helikopter_4_M, 0,
                sizeof(RT_MODEL_helikopter_4_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helikopter_4_M->solverInfo,
                          &helikopter_4_M->Timing.simTimeStep);
    rtsiSetTPtr(&helikopter_4_M->solverInfo, &rtmGetTPtr(helikopter_4_M));
    rtsiSetStepSizePtr(&helikopter_4_M->solverInfo,
                       &helikopter_4_M->Timing.stepSize0);
    rtsiSetdXPtr(&helikopter_4_M->solverInfo, &helikopter_4_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helikopter_4_M->solverInfo, (real_T **)
                         &helikopter_4_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helikopter_4_M->solverInfo,
      &helikopter_4_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helikopter_4_M->solverInfo, (&rtmGetErrorStatus
      (helikopter_4_M)));
    rtsiSetRTModelPtr(&helikopter_4_M->solverInfo, helikopter_4_M);
  }

  rtsiSetSimTimeStep(&helikopter_4_M->solverInfo, MAJOR_TIME_STEP);
  helikopter_4_M->ModelData.intgData.f[0] = helikopter_4_M->ModelData.odeF[0];
  helikopter_4_M->ModelData.contStates = ((real_T *) &helikopter_4_X);
  rtsiSetSolverData(&helikopter_4_M->solverInfo, (void *)
                    &helikopter_4_M->ModelData.intgData);
  rtsiSetSolverName(&helikopter_4_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helikopter_4_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helikopter_4_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helikopter_4_M->Timing.sampleTimes =
      (&helikopter_4_M->Timing.sampleTimesArray[0]);
    helikopter_4_M->Timing.offsetTimes =
      (&helikopter_4_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helikopter_4_M->Timing.sampleTimes[0] = (0.0);
    helikopter_4_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    helikopter_4_M->Timing.offsetTimes[0] = (0.0);
    helikopter_4_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helikopter_4_M, &helikopter_4_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helikopter_4_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helikopter_4_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helikopter_4_M, 35.0);
  helikopter_4_M->Timing.stepSize0 = 0.002;
  helikopter_4_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  helikopter_4_M->Sizes.checksums[0] = (1557197245U);
  helikopter_4_M->Sizes.checksums[1] = (812293119U);
  helikopter_4_M->Sizes.checksums[2] = (72724192U);
  helikopter_4_M->Sizes.checksums[3] = (242791730U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helikopter_4_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helikopter_4_M->extModeInfo,
      &helikopter_4_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helikopter_4_M->extModeInfo,
                        helikopter_4_M->Sizes.checksums);
    rteiSetTPtr(helikopter_4_M->extModeInfo, rtmGetTPtr(helikopter_4_M));
  }

  helikopter_4_M->solverInfoPtr = (&helikopter_4_M->solverInfo);
  helikopter_4_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helikopter_4_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helikopter_4_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helikopter_4_M->ModelData.blockIO = ((void *) &helikopter_4_B);

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_4_B.TmpSignalConversionAtToWorkspac[i] = 0.0;
    }

    helikopter_4_B.u[0] = 0.0;
    helikopter_4_B.u[1] = 0.0;
    helikopter_4_B.TravelCounttorad = 0.0;
    helikopter_4_B.Gain = 0.0;
    helikopter_4_B.travel = 0.0;
    helikopter_4_B.Gain_d = 0.0;
    helikopter_4_B.PitchCounttorad = 0.0;
    helikopter_4_B.Gain_i = 0.0;
    helikopter_4_B.Gain_b = 0.0;
    helikopter_4_B.ElevationCounttorad = 0.0;
    helikopter_4_B.Gain_e = 0.0;
    helikopter_4_B.Sum = 0.0;
    helikopter_4_B.Gain_dg = 0.0;
    helikopter_4_B.Sum1 = 0.0;
    helikopter_4_B.Sum2 = 0.0;
    helikopter_4_B.K_ei = 0.0;
    helikopter_4_B.Gain_l = 0.0;
    helikopter_4_B.BackmotorSaturation = 0.0;
    helikopter_4_B.FrontmotorSaturation = 0.0;
  }

  /* parameters */
  helikopter_4_M->ModelData.defaultParam = ((real_T *)&helikopter_4_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helikopter_4_X;
    helikopter_4_M->ModelData.contStates = (x);
    (void) memset((void *)&helikopter_4_X, 0,
                  sizeof(X_helikopter_4_T));
  }

  /* states (dwork) */
  helikopter_4_M->ModelData.dwork = ((void *) &helikopter_4_DW);
  (void) memset((void *)&helikopter_4_DW, 0,
                sizeof(DW_helikopter_4_T));

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_4_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_4_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_4_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_4_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_4_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_4_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_4_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_4_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helikopter_4_DW.TimeStampA = 0.0;
  helikopter_4_DW.LastUAtTimeA = 0.0;
  helikopter_4_DW.TimeStampB = 0.0;
  helikopter_4_DW.LastUAtTimeB = 0.0;
  helikopter_4_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helikopter_4_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helikopter_4_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  helikopter_4_M->Sizes.numContStates = (4);/* Number of continuous states */
  helikopter_4_M->Sizes.numY = (0);    /* Number of model outputs */
  helikopter_4_M->Sizes.numU = (0);    /* Number of model inputs */
  helikopter_4_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helikopter_4_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helikopter_4_M->Sizes.numBlocks = (57);/* Number of blocks */
  helikopter_4_M->Sizes.numBlockIO = (19);/* Number of block outputs */
  helikopter_4_M->Sizes.numBlockPrms = (141);/* Sum of parameter "widths" */
  return helikopter_4_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
