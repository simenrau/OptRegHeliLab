/*
 * helikopter_3.c
 *
 * Code generation for model "helikopter_3".
 *
 * Model version              : 1.177
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Tue Feb 20 14:47:39 2018
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#include "helikopter_3.h"
#include "helikopter_3_private.h"
#include "helikopter_3_dt.h"

/* Block signals (auto storage) */
B_helikopter_3_T helikopter_3_B;

/* Continuous states */
X_helikopter_3_T helikopter_3_X;

/* Block states (auto storage) */
DW_helikopter_3_T helikopter_3_DW;

/* Real-time model */
RT_MODEL_helikopter_3_T helikopter_3_M_;
RT_MODEL_helikopter_3_T *const helikopter_3_M = &helikopter_3_M_;

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
  helikopter_3_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; i++) {
    *x += h * f0[i];
    x++;
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helikopter_3_output(void)
{
  /* local block i/o variables */
  real_T rtb_Sum3_j[4];
  real_T rtb_Derivative;
  real_T rtb_Backgain;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T *lastU;
  real_T rtb_Gain1_e_idx_4;
  real_T rtb_Gain1_e_idx_5;
  if (rtmIsMajorTimeStep(helikopter_3_M)) {
    /* set solver stop time */
    if (!(helikopter_3_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helikopter_3_M->solverInfo,
                            ((helikopter_3_M->Timing.clockTickH0 + 1) *
        helikopter_3_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helikopter_3_M->solverInfo,
                            ((helikopter_3_M->Timing.clockTick0 + 1) *
        helikopter_3_M->Timing.stepSize0 + helikopter_3_M->Timing.clockTickH0 *
        helikopter_3_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helikopter_3_M)) {
    helikopter_3_M->Timing.t[0] = rtsiGetT(&helikopter_3_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helikopter_3_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S5>/HIL Read Encoder Timebase' */

    /* S-Function Block: helikopter_3/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder(helikopter_3_DW.HILReadEncoderTimebase_Task,
        1, &helikopter_3_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helikopter_3_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helikopter_3_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          helikopter_3_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace' */
  {
    real_T *pDataValues = (real_T *) helikopter_3_DW.FromWorkspace_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) helikopter_3_DW.FromWorkspace_PWORK.TimePtr;
    int_T currTimeIndex = helikopter_3_DW.FromWorkspace_IWORK.PrevIndex;
    real_T t = helikopter_3_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
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

    helikopter_3_DW.FromWorkspace_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          rtb_Derivative = pDataValues[currTimeIndex];
        } else {
          rtb_Derivative = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        rtb_Derivative = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace1' */
  {
    real_T *pDataValues = (real_T *)
      helikopter_3_DW.FromWorkspace1_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helikopter_3_DW.FromWorkspace1_PWORK.TimePtr;
    int_T currTimeIndex = helikopter_3_DW.FromWorkspace1_IWORK.PrevIndex;
    real_T t = helikopter_3_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
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

    helikopter_3_DW.FromWorkspace1_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&rtb_Sum3_j[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 141;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&rtb_Sum3_j[0])[elIdx] = pDataValues[currTimeIndex + 1];
              pDataValues += 141;
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
          for (elIdx = 0; elIdx < 4; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_Sum3_j[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1, f2);
            pDataValues += 141;
          }
        }
      }
    }
  }

  if (rtmIsMajorTimeStep(helikopter_3_M)) {
    /* Gain: '<S5>/Travel: Count to rad' */
    helikopter_3_B.TravelCounttorad = helikopter_3_P.TravelCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o1;

    /* Gain: '<S12>/Gain' */
    helikopter_3_B.Gain = helikopter_3_P.Gain_Gain *
      helikopter_3_B.TravelCounttorad;

    /* Sum: '<Root>/Sum5' incorporates:
     *  Constant: '<Root>/Constant'
     */
    helikopter_3_B.travel = helikopter_3_P.Constant_Value + helikopter_3_B.Gain;
  }

  /* TransferFcn: '<S5>/Travel: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helikopter_3_P.TravelTransferFcn_C *
    helikopter_3_X.TravelTransferFcn_CSTATE;
  rtb_Backgain += helikopter_3_P.TravelTransferFcn_D *
    helikopter_3_B.TravelCounttorad;

  /* Gain: '<S13>/Gain' */
  helikopter_3_B.Gain_d = helikopter_3_P.Gain_Gain_l * rtb_Backgain;
  if (rtmIsMajorTimeStep(helikopter_3_M)) {
    /* Gain: '<S5>/Pitch: Count to rad' */
    helikopter_3_B.PitchCounttorad = helikopter_3_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S9>/Gain' */
    helikopter_3_B.Gain_i = helikopter_3_P.Gain_Gain_a *
      helikopter_3_B.PitchCounttorad;
  }

  /* TransferFcn: '<S5>/Pitch: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helikopter_3_P.PitchTransferFcn_C *
    helikopter_3_X.PitchTransferFcn_CSTATE;
  rtb_Backgain += helikopter_3_P.PitchTransferFcn_D *
    helikopter_3_B.PitchCounttorad;

  /* Gain: '<S10>/Gain' */
  helikopter_3_B.Gain_b = helikopter_3_P.Gain_Gain_ae * rtb_Backgain;

  /* Sum: '<Root>/Sum3' incorporates:
   *  Gain: '<S3>/Gain1'
   */
  rtb_Sum3_j[0] = helikopter_3_P.Gain1_Gain * helikopter_3_B.travel -
    rtb_Sum3_j[0];
  rtb_Sum3_j[1] = helikopter_3_P.Gain1_Gain * helikopter_3_B.Gain_d -
    rtb_Sum3_j[1];
  rtb_Sum3_j[2] = helikopter_3_P.Gain1_Gain * helikopter_3_B.Gain_i -
    rtb_Sum3_j[2];
  rtb_Sum3_j[3] = helikopter_3_P.Gain1_Gain * helikopter_3_B.Gain_b -
    rtb_Sum3_j[3];

  /* Gain: '<Root>/Gain' */
  rtb_Backgain = ((helikopter_3_P.K[0] * rtb_Sum3_j[0] + helikopter_3_P.K[1] *
                   rtb_Sum3_j[1]) + helikopter_3_P.K[2] * rtb_Sum3_j[2]) +
    helikopter_3_P.K[3] * rtb_Sum3_j[3];

  /* Sum: '<Root>/Sum4' */
  helikopter_3_B.u = rtb_Derivative - rtb_Backgain;
  if (rtmIsMajorTimeStep(helikopter_3_M)) {
    /* SignalConversion: '<Root>/TmpSignal ConversionAtTo WorkspaceInport1' */
    helikopter_3_B.TmpSignalConversionAtToWorkspac[0] = helikopter_3_B.u;
    helikopter_3_B.TmpSignalConversionAtToWorkspac[1] = helikopter_3_B.travel;
    helikopter_3_B.TmpSignalConversionAtToWorkspac[2] = helikopter_3_B.Gain_d;
    helikopter_3_B.TmpSignalConversionAtToWorkspac[3] = helikopter_3_B.Gain_i;
    helikopter_3_B.TmpSignalConversionAtToWorkspac[4] = helikopter_3_B.Gain_b;

    /* Gain: '<S5>/Elevation: Count to rad' */
    helikopter_3_B.ElevationCounttorad = helikopter_3_P.ElevationCounttorad_Gain
      * rtb_HILReadEncoderTimebase_o3;

    /* Gain: '<S7>/Gain' */
    helikopter_3_B.Gain_e = helikopter_3_P.Gain_Gain_lv *
      helikopter_3_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helikopter_3_B.Sum = helikopter_3_B.Gain_e +
      helikopter_3_P.elavation_offsetdeg_Value;
  }

  /* TransferFcn: '<S5>/Elevation: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helikopter_3_P.ElevationTransferFcn_C *
    helikopter_3_X.ElevationTransferFcn_CSTATE;
  rtb_Backgain += helikopter_3_P.ElevationTransferFcn_D *
    helikopter_3_B.ElevationCounttorad;

  /* Gain: '<S8>/Gain' */
  helikopter_3_B.Gain_dg = helikopter_3_P.Gain_Gain_n * rtb_Backgain;

  /* Gain: '<S2>/Gain1' */
  rtb_Gain1_e_idx_4 = helikopter_3_P.Gain1_Gain_f * helikopter_3_B.Sum;
  rtb_Gain1_e_idx_5 = helikopter_3_P.Gain1_Gain_f * helikopter_3_B.Gain_dg;

  /* Sum: '<Root>/Sum1' incorporates:
   *  Constant: '<Root>/Vd_bias'
   *  Gain: '<S2>/Gain1'
   *  Gain: '<S6>/K_pd'
   *  Gain: '<S6>/K_pp'
   *  Sum: '<S6>/Sum2'
   *  Sum: '<S6>/Sum3'
   */
  helikopter_3_B.Sum1 = ((helikopter_3_B.u - helikopter_3_P.Gain1_Gain_f *
    helikopter_3_B.Gain_i) * helikopter_3_P.K_pp - helikopter_3_P.Gain1_Gain_f *
    helikopter_3_B.Gain_b * helikopter_3_P.K_pd) + helikopter_3_P.Vd_ff;
  if (rtmIsMajorTimeStep(helikopter_3_M)) {
  }

  /* Integrator: '<S4>/Integrator'
   *
   * Regarding '<S4>/Integrator':
   *  Limited Integrator
   */
  if (helikopter_3_X.Integrator_CSTATE >= helikopter_3_P.Integrator_UpperSat ) {
    helikopter_3_X.Integrator_CSTATE = helikopter_3_P.Integrator_UpperSat;
  } else if (helikopter_3_X.Integrator_CSTATE <=
             (helikopter_3_P.Integrator_LowerSat) ) {
    helikopter_3_X.Integrator_CSTATE = (helikopter_3_P.Integrator_LowerSat);
  }

  rtb_Backgain = helikopter_3_X.Integrator_CSTATE;

  /* Sum: '<S4>/Sum' incorporates:
   *  Constant: '<Root>/elevation_ref'
   */
  rtb_Derivative = helikopter_3_P.elevation_ref_Value - rtb_Gain1_e_idx_4;

  /* Sum: '<Root>/Sum2' incorporates:
   *  Constant: '<Root>/Vs_bias'
   *  Gain: '<S4>/K_ed'
   *  Gain: '<S4>/K_ep'
   *  Sum: '<S4>/Sum1'
   */
  helikopter_3_B.Sum2 = ((helikopter_3_P.K_ep * rtb_Derivative + rtb_Backgain) -
    helikopter_3_P.K_ed * rtb_Gain1_e_idx_5) + helikopter_3_P.Vs_ff;
  if (rtmIsMajorTimeStep(helikopter_3_M)) {
  }

  /* Gain: '<S1>/Back gain' incorporates:
   *  Sum: '<S1>/Subtract'
   */
  rtb_Backgain = (helikopter_3_B.Sum2 - helikopter_3_B.Sum1) *
    helikopter_3_P.Backgain_Gain;

  /* Gain: '<S4>/K_ei' */
  helikopter_3_B.K_ei = helikopter_3_P.K_ei * rtb_Derivative;
  if (rtmIsMajorTimeStep(helikopter_3_M)) {
  }

  /* Derivative: '<S5>/Derivative' */
  if ((helikopter_3_DW.TimeStampA >= helikopter_3_M->Timing.t[0]) &&
      (helikopter_3_DW.TimeStampB >= helikopter_3_M->Timing.t[0])) {
    rtb_Derivative = 0.0;
  } else {
    rtb_Gain1_e_idx_4 = helikopter_3_DW.TimeStampA;
    lastU = &helikopter_3_DW.LastUAtTimeA;
    if (helikopter_3_DW.TimeStampA < helikopter_3_DW.TimeStampB) {
      if (helikopter_3_DW.TimeStampB < helikopter_3_M->Timing.t[0]) {
        rtb_Gain1_e_idx_4 = helikopter_3_DW.TimeStampB;
        lastU = &helikopter_3_DW.LastUAtTimeB;
      }
    } else {
      if (helikopter_3_DW.TimeStampA >= helikopter_3_M->Timing.t[0]) {
        rtb_Gain1_e_idx_4 = helikopter_3_DW.TimeStampB;
        lastU = &helikopter_3_DW.LastUAtTimeB;
      }
    }

    rtb_Derivative = (helikopter_3_B.PitchCounttorad - *lastU) /
      (helikopter_3_M->Timing.t[0] - rtb_Gain1_e_idx_4);
  }

  /* End of Derivative: '<S5>/Derivative' */

  /* Gain: '<S11>/Gain' */
  helikopter_3_B.Gain_l = helikopter_3_P.Gain_Gain_a1 * rtb_Derivative;
  if (rtmIsMajorTimeStep(helikopter_3_M)) {
  }

  /* Saturate: '<S5>/Back motor: Saturation' */
  if (rtb_Backgain > helikopter_3_P.BackmotorSaturation_UpperSat) {
    helikopter_3_B.BackmotorSaturation =
      helikopter_3_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Backgain < helikopter_3_P.BackmotorSaturation_LowerSat) {
    helikopter_3_B.BackmotorSaturation =
      helikopter_3_P.BackmotorSaturation_LowerSat;
  } else {
    helikopter_3_B.BackmotorSaturation = rtb_Backgain;
  }

  /* End of Saturate: '<S5>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helikopter_3_M)) {
  }

  /* Gain: '<S1>/Front gain' incorporates:
   *  Sum: '<S1>/Add'
   */
  rtb_Gain1_e_idx_4 = (helikopter_3_B.Sum1 + helikopter_3_B.Sum2) *
    helikopter_3_P.Frontgain_Gain;

  /* Saturate: '<S5>/Front motor: Saturation' */
  if (rtb_Gain1_e_idx_4 > helikopter_3_P.FrontmotorSaturation_UpperSat) {
    helikopter_3_B.FrontmotorSaturation =
      helikopter_3_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Gain1_e_idx_4 < helikopter_3_P.FrontmotorSaturation_LowerSat) {
    helikopter_3_B.FrontmotorSaturation =
      helikopter_3_P.FrontmotorSaturation_LowerSat;
  } else {
    helikopter_3_B.FrontmotorSaturation = rtb_Gain1_e_idx_4;
  }

  /* End of Saturate: '<S5>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helikopter_3_M)) {
    /* S-Function (hil_write_analog_block): '<S5>/HIL Write Analog' */

    /* S-Function Block: helikopter_3/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helikopter_3_DW.HILWriteAnalog_Buffer[0] =
        helikopter_3_B.FrontmotorSaturation;
      helikopter_3_DW.HILWriteAnalog_Buffer[1] =
        helikopter_3_B.BackmotorSaturation;
      result = hil_write_analog(helikopter_3_DW.HILInitialize_Card,
        helikopter_3_P.HILWriteAnalog_channels, 2,
        &helikopter_3_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void helikopter_3_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S5>/Derivative' */
  if (helikopter_3_DW.TimeStampA == (rtInf)) {
    helikopter_3_DW.TimeStampA = helikopter_3_M->Timing.t[0];
    lastU = &helikopter_3_DW.LastUAtTimeA;
  } else if (helikopter_3_DW.TimeStampB == (rtInf)) {
    helikopter_3_DW.TimeStampB = helikopter_3_M->Timing.t[0];
    lastU = &helikopter_3_DW.LastUAtTimeB;
  } else if (helikopter_3_DW.TimeStampA < helikopter_3_DW.TimeStampB) {
    helikopter_3_DW.TimeStampA = helikopter_3_M->Timing.t[0];
    lastU = &helikopter_3_DW.LastUAtTimeA;
  } else {
    helikopter_3_DW.TimeStampB = helikopter_3_M->Timing.t[0];
    lastU = &helikopter_3_DW.LastUAtTimeB;
  }

  *lastU = helikopter_3_B.PitchCounttorad;

  /* End of Update for Derivative: '<S5>/Derivative' */
  if (rtmIsMajorTimeStep(helikopter_3_M)) {
    rt_ertODEUpdateContinuousStates(&helikopter_3_M->solverInfo);
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
  if (!(++helikopter_3_M->Timing.clockTick0)) {
    ++helikopter_3_M->Timing.clockTickH0;
  }

  helikopter_3_M->Timing.t[0] = rtsiGetSolverStopTime
    (&helikopter_3_M->solverInfo);

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
    if (!(++helikopter_3_M->Timing.clockTick1)) {
      ++helikopter_3_M->Timing.clockTickH1;
    }

    helikopter_3_M->Timing.t[1] = helikopter_3_M->Timing.clockTick1 *
      helikopter_3_M->Timing.stepSize1 + helikopter_3_M->Timing.clockTickH1 *
      helikopter_3_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void helikopter_3_derivatives(void)
{
  XDot_helikopter_3_T *_rtXdot;
  _rtXdot = ((XDot_helikopter_3_T *) helikopter_3_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S5>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += helikopter_3_P.TravelTransferFcn_A *
    helikopter_3_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helikopter_3_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S5>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += helikopter_3_P.PitchTransferFcn_A *
    helikopter_3_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helikopter_3_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S5>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE += helikopter_3_P.ElevationTransferFcn_A *
    helikopter_3_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += helikopter_3_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S4>/Integrator' */
  {
    boolean_T lsat;
    boolean_T usat;
    lsat = ( helikopter_3_X.Integrator_CSTATE <=
            (helikopter_3_P.Integrator_LowerSat) );
    usat = ( helikopter_3_X.Integrator_CSTATE >=
            helikopter_3_P.Integrator_UpperSat );
    if ((!lsat && !usat) ||
        (lsat && (helikopter_3_B.K_ei > 0)) ||
        (usat && (helikopter_3_B.K_ei < 0)) ) {
      ((XDot_helikopter_3_T *) helikopter_3_M->ModelData.derivs)
        ->Integrator_CSTATE = helikopter_3_B.K_ei;
    } else {
      /* in saturation */
      ((XDot_helikopter_3_T *) helikopter_3_M->ModelData.derivs)
        ->Integrator_CSTATE = 0.0;
    }
  }
}

/* Model initialize function */
void helikopter_3_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopter_3/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &helikopter_3_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options(helikopter_3_DW.HILInitialize_Card,
      "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helikopter_3_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
      return;
    }

    if ((helikopter_3_P.HILInitialize_set_analog_input_ && !is_switching) ||
        (helikopter_3_P.HILInitialize_set_analog_inpu_m && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helikopter_3_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = helikopter_3_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helikopter_3_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = helikopter_3_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges(helikopter_3_DW.HILInitialize_Card,
        helikopter_3_P.HILInitialize_analog_input_chan, 8U,
        &helikopter_3_DW.HILInitialize_AIMinimums[0],
        &helikopter_3_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter_3_P.HILInitialize_set_analog_output && !is_switching) ||
        (helikopter_3_P.HILInitialize_set_analog_outp_b && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &helikopter_3_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = helikopter_3_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &helikopter_3_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = helikopter_3_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges(helikopter_3_DW.HILInitialize_Card,
        helikopter_3_P.HILInitialize_analog_output_cha, 8U,
        &helikopter_3_DW.HILInitialize_AOMinimums[0],
        &helikopter_3_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter_3_P.HILInitialize_set_analog_outp_e && !is_switching) ||
        (helikopter_3_P.HILInitialize_set_analog_outp_j && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helikopter_3_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helikopter_3_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(helikopter_3_DW.HILInitialize_Card,
        helikopter_3_P.HILInitialize_analog_output_cha, 8U,
        &helikopter_3_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
        return;
      }
    }

    if (helikopter_3_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helikopter_3_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helikopter_3_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helikopter_3_DW.HILInitialize_Card,
         helikopter_3_P.HILInitialize_analog_output_cha, 8U,
         &helikopter_3_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter_3_P.HILInitialize_set_encoder_param && !is_switching) ||
        (helikopter_3_P.HILInitialize_set_encoder_par_m && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helikopter_3_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = helikopter_3_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode
        (helikopter_3_DW.HILInitialize_Card,
         helikopter_3_P.HILInitialize_encoder_channels, 8U,
         (t_encoder_quadrature_mode *)
         &helikopter_3_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter_3_P.HILInitialize_set_encoder_count && !is_switching) ||
        (helikopter_3_P.HILInitialize_set_encoder_cou_k && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helikopter_3_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] =
            helikopter_3_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts(helikopter_3_DW.HILInitialize_Card,
        helikopter_3_P.HILInitialize_encoder_channels, 8U,
        &helikopter_3_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter_3_P.HILInitialize_set_pwm_params_at && !is_switching) ||
        (helikopter_3_P.HILInitialize_set_pwm_params__f && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues = &helikopter_3_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helikopter_3_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(helikopter_3_DW.HILInitialize_Card,
        helikopter_3_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &helikopter_3_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          helikopter_3_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues = &helikopter_3_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            helikopter_3_DW.HILInitialize_POSortedChans[num_duty_cycle_modes] =
              p_HILInitialize_pwm_channels[i1];
            helikopter_3_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes] =
              helikopter_3_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            helikopter_3_DW.HILInitialize_POSortedChans[7U - num_frequency_modes]
              = p_HILInitialize_pwm_channels[i1];
            helikopter_3_DW.HILInitialize_POSortedFreqs[7U - num_frequency_modes]
              = helikopter_3_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(helikopter_3_DW.HILInitialize_Card,
          &helikopter_3_DW.HILInitialize_POSortedChans[0], num_duty_cycle_modes,
          &helikopter_3_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(helikopter_3_DW.HILInitialize_Card,
          &helikopter_3_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &helikopter_3_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues = &helikopter_3_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helikopter_3_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues =
          &helikopter_3_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = helikopter_3_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helikopter_3_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = helikopter_3_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration(helikopter_3_DW.HILInitialize_Card,
        helikopter_3_P.HILInitialize_pwm_channels, 8U,
        (t_pwm_configuration *) &helikopter_3_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &helikopter_3_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &helikopter_3_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs = &helikopter_3_DW.HILInitialize_POSortedFreqs
          [0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = helikopter_3_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helikopter_3_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helikopter_3_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband(helikopter_3_DW.HILInitialize_Card,
        helikopter_3_P.HILInitialize_pwm_channels, 8U,
        &helikopter_3_DW.HILInitialize_POSortedFreqs[0],
        &helikopter_3_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
        return;
      }
    }

    if ((helikopter_3_P.HILInitialize_set_pwm_outputs_a && !is_switching) ||
        (helikopter_3_P.HILInitialize_set_pwm_outputs_g && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helikopter_3_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helikopter_3_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(helikopter_3_DW.HILInitialize_Card,
        helikopter_3_P.HILInitialize_pwm_channels, 8U,
        &helikopter_3_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
        return;
      }
    }

    if (helikopter_3_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues = &helikopter_3_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helikopter_3_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helikopter_3_DW.HILInitialize_Card,
         helikopter_3_P.HILInitialize_pwm_channels, 8U,
         &helikopter_3_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S5>/HIL Read Encoder Timebase' */

  /* S-Function Block: helikopter_3/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader(helikopter_3_DW.HILInitialize_Card,
      helikopter_3_P.HILReadEncoderTimebase_samples_,
      helikopter_3_P.HILReadEncoderTimebase_channels, 3,
      &helikopter_3_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
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
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877572486291,
      0.52359877557845225, 0.52359877555326229, 0.52359877555746936,
      0.52359877556562584, 0.52359877557288137, 0.52359877554548584,
      0.52359877553665346, 0.52359877555443346, 0.5235987755206154,
      0.52359877562053325, 0.52359877524968768, 0.52359877484315143,
      0.040518810698263567, -0.52359877437196545, -0.52359877512717179,
      -0.52359877526592147, -0.5235987753134721, -0.52359877533487,
      -0.523598775345218, -0.52359877535044252, -0.52359877537617627,
      -0.52359877548306333, -0.523598775506603, -0.52359877489184214,
      -0.52359877481970152, -0.52359877557068168, -0.52359877094471263,
      -0.52359877551985834, -0.52359877539517929, -0.27115064354188412,
      0.25927927578549631, 0.48087695413606063, 0.4987511928136914,
      0.40561121215652673, 0.26963539787265067, 0.13729579402594153,
      0.033346396657273507, -0.033848402410424847, -0.066701083952162413,
      -0.073262755478817221, -0.063255253574726458, -0.045492531082207861,
      -0.026575489037097867, -0.010562581974817102, 0.00073584232019924,
      0.0071591766858969553, 0.0095360033330118225, 0.0091129270543912066,
      0.0071436079979592488, 0.0046486791079056126, 0.0023203675973239571,
      0.00052921657033308664, -0.000611160907639516, -0.0011578707070388685,
      -0.0012563722288431124, -0.0010759271319775152, -0.00076760520632160182,
      -0.0004433351102883235, -0.00017114074075662592, 1.932068865892592E-5,
      0.0001262849020432205, 0.00016452135952978387, 0.00015560844490730902,
      0.00012099740589228226, 7.8007986287560263E-5, 3.8284219858993053E-5,
      7.9711677931131716E-6, -1.1142262079942415E-5, -2.0135574553800214E-5,
      -2.1550035145341683E-5, -1.829519256963999E-5, -1.2944130129208504E-5,
      -7.3877966160931942E-6, -2.7632739406515405E-6, 4.452860843256655E-7,
      2.2245329538871646E-6, 2.8370838953546942E-6, 2.6561652860161438E-6,
      2.0486508611641125E-6, 1.3083311124929351E-6, 6.3097706193751172E-7,
      1.1839227907062391E-7, -2.0146773020816847E-7, -3.4885130660681846E-7,
      -3.682074479057346E-7, -3.0967161878601231E-7, -2.17232691212071E-7,
      -1.2300503819499287E-7, -4.6113297706103258E-8, 5.5489044264976952E-9,
      3.2232480262311778E-8, 3.90511143881656E-8, 3.3171028432197061E-8,
      2.1733080440470294E-8, 1.04823846481439E-8, 2.9427562064199426E-9,
      1.0910346494880082E-12, 1.17420506460252E-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0 } ;

    helikopter_3_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    helikopter_3_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    helikopter_3_DW.FromWorkspace_IWORK.PrevIndex = 0;
  }

  /* Start for FromWorkspace: '<Root>/From Workspace1' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1378421413117508, 3.1262155534054883,
      3.1033092999754026, 3.0666274151347692, 3.0144539223364673,
      2.9456562770593426, 2.8595077632358281, 2.7555515879093129,
      2.633505110437834, 2.4931956059848517, 2.3345185760236493,
      2.1574113214405486, 1.9618364840177469, 1.7512322141750956,
      1.5334325246981109, 1.3160147891542651, 1.1049443296011534,
      0.90449084479831, 0.7175464934589536, 0.54600106629293976,
      0.39105845507039849, 0.25347127904892558, 0.13370324215217794,
      0.032036751370923877, -0.0513578160463011, -0.1163784732413321,
      -0.16296457378249049, -0.19108025765266029, -0.20070442434582511,
      -0.19182471096143758, -0.16624217474007813, -0.12974116179696832,
      -0.089726926293523551, -0.052723515398272289, -0.023000134521151121,
      -0.0023053972876263111, 0.0096765980116060961, 0.01456688328932939,
      0.014534346416980504, 0.011709724038464348, 0.0078277388107604769,
      0.0040803545787843074, 0.0011267444623190455, -0.000804072009919818,
      -0.0017742728892691947, -0.0020024624960880131, -0.001756306749768666,
      -0.0012811354468802185, -0.000762811948497394, -0.00031746755298479596,
      1.2480308738537328E-6, 0.00018612742178554496, 0.00025830146377552217,
      0.0002513627857909213, 0.00019978725470782492, 0.00013204211538005766,
      6.7719868080089534E-5, 1.7552426779631045E-5, -1.4904565203088609E-5,
      -3.0936019049990863E-5, -3.4395184998137732E-5, -2.9898132065164876E-5,
      -2.1630813834488025E-5, -1.2737575601807046E-5, -5.1630290665359674E-6,
      2.1284844807861584E-7, 3.29486630692465E-6, 4.4616388498166359E-6,
      4.2953455739648487E-6, 3.3862717425768869E-6, 2.2177535587892247E-6,
      1.1195852769568215E-6, 2.7000197885472171E-7, -2.7450025065164281E-7,
      -5.3881705848029433E-7, -5.9040647925762665E-7, -5.0868607552740728E-7,
      -3.6499678282122473E-7, -2.1249390647110511E-7, -8.3716053456366212E-8,
      6.9266622863088569E-9, 5.827789772271476E-8, 7.7096354521161139E-8,
      7.3457539768767517E-8, 5.7455705301552647E-8, 3.7285334331559308E-8,
      1.8487557952232465E-8, 4.0166537117901652E-9, -5.2288132776261558E-9,
      -9.7048620989825287E-9, -1.0554983497986504E-8, -9.0902046610952487E-9,
      -6.4555910830642E-9, -3.4720103182822809E-9, -6.1410490659793015E-10,
      1.9249189728973376E-9, 4.1301198784909E-9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.015002048889369925,
      -0.046506351601760525, -0.0916250136970537, -0.14672753933924165,
      -0.20869397116991958, -0.27519058108520866, -0.34459405527076686,
      -0.41582470128277321, -0.48818590986262544, -0.5612380177886388,
      -0.63470811982151876, -0.70842901830911265, -0.78229934966791714,
      -0.84241707934731547, -0.87119875788464862, -0.8696709421520924,
      -0.84428183818915714, -0.801813939188084, -0.74777740533413584,
      -0.68618170864076544, -0.61977044486687538, -0.55034870406260172,
      -0.47907214756370081, -0.40666596310172642, -0.33357826964561005,
      -0.26008262875683419, -0.18634440214134379, -0.11246273545738936,
      -0.03849666674936944, 0.0355188535608399, 0.10233014490872774,
      0.14600405179572912, 0.16005694203706888, 0.14801364360429492,
      0.11889352353177451, 0.08277894895738909, 0.047927981220219476,
      0.01956114113418302, -0.000130147466105697, -0.011298489490774775,
      -0.015527940887525639, -0.014989536904614827, -0.011814440442571203,
      -0.0077232658656656078, -0.0038808034941076605, -0.00091275840398542839,
      0.00098462300856723417, 0.001900685234843637, 0.0020732940168211437,
      0.0017813776053402382, 0.0012748623587244448, 0.00073951758693661094,
      0.00028869619124975487, -2.7754688648557442E-5, -0.00020630210104253952,
      -0.00027098053402122306, -0.00025728896591002649, -0.00020066974191198793,
      -0.00012982794464103263, -6.4125792097763018E-5, -1.3836640502741463E-5,
      1.7988235021737424E-5, 3.3069296212553412E-5, 3.5572976220569921E-5,
      3.0298209430930314E-5, 2.1503533348304338E-5, 1.232809472523014E-5,
      4.6671134614139469E-6, -6.6514981356114446E-7, -3.636272035705845E-6,
      -4.6740494453046457E-6, -4.3926498374836079E-6, -3.3983099025623958E-6,
      -2.1779856281794546E-6, -1.0572439414686026E-6, -2.0633439326332539E-7,
      3.2690490476688129E-7, 5.7478046067073378E-7, 6.1003479524648213E-7,
      5.1513470190495924E-7, 3.6259415281670392E-7, 2.0542823159162726E-7,
      7.5297117039789178E-8, -1.4531969163570834E-8, -6.3984048022855849E-8,
      -8.06581940339697E-8, -7.5167815671303726E-8, -5.7860327115765553E-8,
      -3.6958578111661625E-8, -1.7880905439421859E-8, -3.377195750012252E-9,
      5.882405193568657E-9, 1.0561744158127849E-8, 1.195761290513133E-8,
      1.1454911492741042E-8, 1.0179385363984717E-8, 8.8440934683778891E-9,
      7.7609501756476934E-9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.10602875208525725, 0.22266037934805066, 0.31888147182660104,
      0.38944360630686331, 0.43795507376204207, 0.46997264228392854,
      0.4905172487472646, 0.503431001377275, 0.51142138581872776,
      0.51630439852689947, 0.51925862123258226, 0.52103115478951523,
      0.52208728925134873, 0.42488915301905361, 0.20341791151210295,
      -0.01079801808820275, -0.17944049012078328, -0.30014689066215833,
      -0.38190958347214837, -0.435334859401467, -0.46936944834975347,
      -0.49064635022946074, -0.50375547916192676, -0.51173925817379462,
      -0.51655590346028935, -0.51943912010973492, -0.521153623377034,
      -0.52216740299784681, -0.52276392426885, -0.52311342930641025,
      -0.47219669046295565, -0.30867049379749911, -0.099320461099652632,
      0.085117433708898113, 0.2058098870241295, 0.25524402009881253,
      0.24631332957414778, 0.20048599177476964, 0.13917050726795585,
      0.078933575982570725, 0.029892147147810667, -0.0038052337210180182,
      -0.022440369149479639, -0.028914859393831972, -0.027157056515398043,
      -0.020977009130463325, -0.013409967168707619, -0.0064743779486783294,
      -0.0012199329464978885, 0.0020631537079142653, 0.0035798563143721895,
      0.0037836123876234442, 0.0031862334467806517, 0.0022365539604056223,
      0.0012619049206196665, 0.00045712246331820224, -9.6766773062873469E-5,
      -0.00040016304980304139, -0.00050068276576521178, -0.00046435771977952969,
      -0.00035542451566007363, -0.00022492606432902111, -0.0001065871794699697,
      -1.7695052863707707E-5, 3.7280038332774511E-5, 6.21574128805545E-5,
      6.4848497094391745E-5, 5.4144890885568606E-5, 3.7686401531883621E-5,
      2.0998758345042938E-5, 7.3346155992207474E-6, -1.988823735584765E-6,
      -7.0276136790722405E-6, -8.6247846817239445E-6, -7.9209729849589229E-6,
      -6.0139024851573751E-6, -3.7687304825608873E-6, -1.7518885402475444E-6,
      -2.4916295287936134E-7, 6.7071883464331741E-7, 1.0780995616334753E-6,
      1.1107899253595653E-6, 9.1971823701319692E-7, 6.3487889292609408E-7,
      3.4950962593433063E-7, 1.1784771521670897E-7, -3.880267791092265E-8,
      -1.2232149952551376E-7, -1.4772435964449708E-7, -1.348324847879648E-7,
      -1.0250546132793513E-7, -6.5442103559251149E-8, -3.3070561199872354E-8,
      -9.8642477356662142E-9, 3.554121042923788E-9, 9.0161509785043366E-9,
      9.4385525926703615E-9, 7.6552431154563572E-9, 5.5658430442447921E-9, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.42411500835944105,
      0.46652650907446352, 0.38488436993749131, 0.28224853794433913,
      0.19404586984400488, 0.12807027411083594, 0.082178425876633973,
      0.051655010543331416, 0.031961537789100811, 0.01953205085597667,
      0.011816890846020965, 0.0070901342510213926, 0.0042245378706239935,
      -0.3887925449058906, -0.88588496600451272, -0.856863718377933,
      -0.67456988810703222, -0.48282560214221054, -0.32705077121667042,
      -0.2137011036939846, -0.13613835576985586, -0.08510760749553925,
      -0.052436515706574066, -0.031935116024181788, -0.019266581122688881,
      -0.011532866574492531, -0.0068580130459062289, -0.0040551184599617162,
      -0.0023860850607230716, -0.0013980201269508292, 0.20366695539710813,
      0.654104786685116, 0.83740013081467568, 0.7377515792574928,
      0.48276981328421542, 0.19773653232202199, -0.0357227620753692,
      -0.1833093511742227, -0.24526193800396523, -0.24094772511825063,
      -0.19616571531575042, -0.13478952345202491, -0.074540541690556641,
      -0.025897960954119479, 0.0070312115370255671, 0.0247201895630287,
      0.030268167870312682, 0.027742356903407005, 0.021017780032011612,
      0.013132346640938462, 0.0060668104491215426, 0.00081502431629486565,
      -0.0023895157400813257, -0.0037987179222102692, -0.0038985961358539768,
      -0.0032191298059160113, -0.0022155569222344568, -0.0012135850836708256,
      -0.00040207884055883569, 0.00014530020723257434, 0.00043573283976767037,
      0.00052199382861405617, 0.00047335556272605164, 0.000355568529714894,
      0.00021990038807577489, 9.9509521480965952E-5, 1.076436014519501E-5,
      -4.2814401545446572E-5, -6.5833934124893924E-5, -6.6750549457516726E-5,
      -5.4656547693442763E-5, -3.7293734049376049E-5, -2.0155136484103895E-5,
      -6.3886607207608136E-6, 2.8152700769060888E-6, 7.628305289052197E-6,
      8.9807113002319551E-6, 8.0673910590993747E-6, 6.0109256393187358E-6,
      3.6795504399367191E-6, 1.6295461978066352E-6, 1.3078474475036416E-7,
      -7.6426346353946984E-7, -1.1393340865024081E-6, -1.14145377812105E-6,
      -9.26624353024483E-7, -6.2657828266452282E-7, -3.3405199661236077E-7,
      -1.0158815062992968E-7, 5.159078927213283E-8, 1.2933138368612236E-7,
      1.4827672092073955E-7, 1.2950945928351881E-7, 9.2848543702828218E-8,
      5.3696764960363651E-8, 2.1871409588325837E-8, 1.7128963026677487E-9,
      -7.1099480628523713E-9, -8.3343104388426182E-9, -7.0341779674039115E-9,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helikopter_3_DW.FromWorkspace1_PWORK.TimePtr = (void *) pTimeValues0;
    helikopter_3_DW.FromWorkspace1_PWORK.DataPtr = (void *) pDataValues0;
    helikopter_3_DW.FromWorkspace1_IWORK.PrevIndex = 0;
  }

  /* InitializeConditions for TransferFcn: '<S5>/Travel: Transfer Fcn' */
  helikopter_3_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S5>/Pitch: Transfer Fcn' */
  helikopter_3_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S5>/Elevation: Transfer Fcn' */
  helikopter_3_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S4>/Integrator' */
  helikopter_3_X.Integrator_CSTATE = helikopter_3_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S5>/Derivative' */
  helikopter_3_DW.TimeStampA = (rtInf);
  helikopter_3_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void helikopter_3_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helikopter_3/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helikopter_3_DW.HILInitialize_Card);
    hil_monitor_stop_all(helikopter_3_DW.HILInitialize_Card);
    is_switching = false;
    if ((helikopter_3_P.HILInitialize_set_analog_out_ex && !is_switching) ||
        (helikopter_3_P.HILInitialize_set_analog_outp_c && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helikopter_3_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helikopter_3_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((helikopter_3_P.HILInitialize_set_pwm_output_ap && !is_switching) ||
        (helikopter_3_P.HILInitialize_set_pwm_outputs_p && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helikopter_3_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helikopter_3_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helikopter_3_DW.HILInitialize_Card
                         , helikopter_3_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , helikopter_3_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helikopter_3_DW.HILInitialize_AOVoltages[0]
                         , &helikopter_3_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(helikopter_3_DW.HILInitialize_Card,
            helikopter_3_P.HILInitialize_analog_output_cha,
            num_final_analog_outputs, &helikopter_3_DW.HILInitialize_AOVoltages
            [0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(helikopter_3_DW.HILInitialize_Card,
            helikopter_3_P.HILInitialize_pwm_channels, num_final_pwm_outputs,
            &helikopter_3_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helikopter_3_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helikopter_3_DW.HILInitialize_Card);
    hil_monitor_delete_all(helikopter_3_DW.HILInitialize_Card);
    hil_close(helikopter_3_DW.HILInitialize_Card);
    helikopter_3_DW.HILInitialize_Card = NULL;
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
  helikopter_3_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  helikopter_3_update();
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
  helikopter_3_initialize();
}

void MdlTerminate(void)
{
  helikopter_3_terminate();
}

/* Registration function */
RT_MODEL_helikopter_3_T *helikopter_3(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helikopter_3_P.Integrator_UpperSat = rtInf;
  helikopter_3_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helikopter_3_M, 0,
                sizeof(RT_MODEL_helikopter_3_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helikopter_3_M->solverInfo,
                          &helikopter_3_M->Timing.simTimeStep);
    rtsiSetTPtr(&helikopter_3_M->solverInfo, &rtmGetTPtr(helikopter_3_M));
    rtsiSetStepSizePtr(&helikopter_3_M->solverInfo,
                       &helikopter_3_M->Timing.stepSize0);
    rtsiSetdXPtr(&helikopter_3_M->solverInfo, &helikopter_3_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helikopter_3_M->solverInfo, (real_T **)
                         &helikopter_3_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helikopter_3_M->solverInfo,
      &helikopter_3_M->Sizes.numContStates);
    rtsiSetErrorStatusPtr(&helikopter_3_M->solverInfo, (&rtmGetErrorStatus
      (helikopter_3_M)));
    rtsiSetRTModelPtr(&helikopter_3_M->solverInfo, helikopter_3_M);
  }

  rtsiSetSimTimeStep(&helikopter_3_M->solverInfo, MAJOR_TIME_STEP);
  helikopter_3_M->ModelData.intgData.f[0] = helikopter_3_M->ModelData.odeF[0];
  helikopter_3_M->ModelData.contStates = ((real_T *) &helikopter_3_X);
  rtsiSetSolverData(&helikopter_3_M->solverInfo, (void *)
                    &helikopter_3_M->ModelData.intgData);
  rtsiSetSolverName(&helikopter_3_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helikopter_3_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helikopter_3_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helikopter_3_M->Timing.sampleTimes =
      (&helikopter_3_M->Timing.sampleTimesArray[0]);
    helikopter_3_M->Timing.offsetTimes =
      (&helikopter_3_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helikopter_3_M->Timing.sampleTimes[0] = (0.0);
    helikopter_3_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    helikopter_3_M->Timing.offsetTimes[0] = (0.0);
    helikopter_3_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helikopter_3_M, &helikopter_3_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helikopter_3_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helikopter_3_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helikopter_3_M, 35.0);
  helikopter_3_M->Timing.stepSize0 = 0.002;
  helikopter_3_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  helikopter_3_M->Sizes.checksums[0] = (490482735U);
  helikopter_3_M->Sizes.checksums[1] = (2593370399U);
  helikopter_3_M->Sizes.checksums[2] = (4049143768U);
  helikopter_3_M->Sizes.checksums[3] = (1355316700U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helikopter_3_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helikopter_3_M->extModeInfo,
      &helikopter_3_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helikopter_3_M->extModeInfo,
                        helikopter_3_M->Sizes.checksums);
    rteiSetTPtr(helikopter_3_M->extModeInfo, rtmGetTPtr(helikopter_3_M));
  }

  helikopter_3_M->solverInfoPtr = (&helikopter_3_M->solverInfo);
  helikopter_3_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helikopter_3_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helikopter_3_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helikopter_3_M->ModelData.blockIO = ((void *) &helikopter_3_B);

  {
    int_T i;
    for (i = 0; i < 5; i++) {
      helikopter_3_B.TmpSignalConversionAtToWorkspac[i] = 0.0;
    }

    helikopter_3_B.TravelCounttorad = 0.0;
    helikopter_3_B.Gain = 0.0;
    helikopter_3_B.travel = 0.0;
    helikopter_3_B.Gain_d = 0.0;
    helikopter_3_B.PitchCounttorad = 0.0;
    helikopter_3_B.Gain_i = 0.0;
    helikopter_3_B.Gain_b = 0.0;
    helikopter_3_B.u = 0.0;
    helikopter_3_B.ElevationCounttorad = 0.0;
    helikopter_3_B.Gain_e = 0.0;
    helikopter_3_B.Sum = 0.0;
    helikopter_3_B.Gain_dg = 0.0;
    helikopter_3_B.Sum1 = 0.0;
    helikopter_3_B.Sum2 = 0.0;
    helikopter_3_B.K_ei = 0.0;
    helikopter_3_B.Gain_l = 0.0;
    helikopter_3_B.BackmotorSaturation = 0.0;
    helikopter_3_B.FrontmotorSaturation = 0.0;
  }

  /* parameters */
  helikopter_3_M->ModelData.defaultParam = ((real_T *)&helikopter_3_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helikopter_3_X;
    helikopter_3_M->ModelData.contStates = (x);
    (void) memset((void *)&helikopter_3_X, 0,
                  sizeof(X_helikopter_3_T));
  }

  /* states (dwork) */
  helikopter_3_M->ModelData.dwork = ((void *) &helikopter_3_DW);
  (void) memset((void *)&helikopter_3_DW, 0,
                sizeof(DW_helikopter_3_T));

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_3_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_3_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_3_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_3_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_3_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_3_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_3_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 8; i++) {
      helikopter_3_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helikopter_3_DW.TimeStampA = 0.0;
  helikopter_3_DW.LastUAtTimeA = 0.0;
  helikopter_3_DW.TimeStampB = 0.0;
  helikopter_3_DW.LastUAtTimeB = 0.0;
  helikopter_3_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helikopter_3_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helikopter_3_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  helikopter_3_M->Sizes.numContStates = (4);/* Number of continuous states */
  helikopter_3_M->Sizes.numY = (0);    /* Number of model outputs */
  helikopter_3_M->Sizes.numU = (0);    /* Number of model inputs */
  helikopter_3_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helikopter_3_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helikopter_3_M->Sizes.numBlocks = (62);/* Number of blocks */
  helikopter_3_M->Sizes.numBlockIO = (19);/* Number of block outputs */
  helikopter_3_M->Sizes.numBlockPrms = (147);/* Sum of parameter "widths" */
  return helikopter_3_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
