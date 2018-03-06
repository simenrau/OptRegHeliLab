/*
 * helikopter_4.c
 *
 * Code generation for model "helikopter_4".
 *
 * Model version              : 1.185
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Tue Mar 06 17:40:03 2018
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
  real_T rtb_Sum3_j[6];
  real_T rtb_Backgain;
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T *lastU;
  real_T tmp[6];
  int32_T i;
  int32_T i_0;
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
              pDataValues += 141;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 2; ++elIdx) {
              (&rtb_FromWorkspace[0])[elIdx] = pDataValues[currTimeIndex + 1];
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
          for (elIdx = 0; elIdx < 2; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_FromWorkspace[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1,
              f2);
            pDataValues += 141;
          }
        }
      }
    }
  }

  /* FromWorkspace: '<Root>/From Workspace1' */
  {
    real_T *pDataValues = (real_T *)
      helikopter_4_DW.FromWorkspace1_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *)
      helikopter_4_DW.FromWorkspace1_PWORK.TimePtr;
    int_T currTimeIndex = helikopter_4_DW.FromWorkspace1_IWORK.PrevIndex;
    real_T t = helikopter_4_M->Timing.t[0];

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

    helikopter_4_DW.FromWorkspace1_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&rtb_Sum3_j[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 141;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
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
          for (elIdx = 0; elIdx < 6; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_Sum3_j[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1, f2);
            pDataValues += 141;
          }
        }
      }
    }
  }

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

    /* Sum: '<Root>/Sum6' incorporates:
     *  Constant: '<Root>/Constant1'
     */
    helikopter_4_B.elevation = helikopter_4_B.Gain_e +
      helikopter_4_P.Constant1_Value;
  }

  /* TransferFcn: '<S5>/Elevation: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helikopter_4_P.ElevationTransferFcn_C *
    helikopter_4_X.ElevationTransferFcn_CSTATE;
  rtb_Backgain += helikopter_4_P.ElevationTransferFcn_D *
    helikopter_4_B.ElevationCounttorad;

  /* Gain: '<S8>/Gain' */
  helikopter_4_B.Gain_dg = helikopter_4_P.Gain_Gain_n * rtb_Backgain;

  /* Gain: '<S3>/Gain1' */
  tmp[0] = helikopter_4_B.travel;
  tmp[1] = helikopter_4_B.Gain_d;
  tmp[2] = helikopter_4_B.Gain_i;
  tmp[3] = helikopter_4_B.Gain_b;
  tmp[4] = helikopter_4_B.elevation;
  tmp[5] = helikopter_4_B.Gain_dg;

  /* Sum: '<Root>/Sum3' incorporates:
   *  Gain: '<S3>/Gain1'
   */
  for (i = 0; i < 6; i++) {
    rtb_Sum3_j[i] = helikopter_4_P.Gain1_Gain * tmp[i] - rtb_Sum3_j[i];
  }

  /* End of Sum: '<Root>/Sum3' */

  /* Sum: '<Root>/Sum4' incorporates:
   *  Gain: '<Root>/Gain'
   */
  for (i = 0; i < 2; i++) {
    rtb_Gain1_idx_4 = 0.0;
    for (i_0 = 0; i_0 < 6; i_0++) {
      rtb_Gain1_idx_4 += helikopter_4_P.K[(i_0 << 1) + i] * rtb_Sum3_j[i_0];
    }

    helikopter_4_B.u[i] = rtb_FromWorkspace[i] - rtb_Gain1_idx_4;
  }

  /* End of Sum: '<Root>/Sum4' */
  if (rtmIsMajorTimeStep(helikopter_4_M)) {
    /* SignalConversion: '<Root>/TmpSignal ConversionAtTo WorkspaceInport1' */
    helikopter_4_B.TmpSignalConversionAtToWorkspac[0] = helikopter_4_B.u[0];
    helikopter_4_B.TmpSignalConversionAtToWorkspac[1] = helikopter_4_B.u[1];
    helikopter_4_B.TmpSignalConversionAtToWorkspac[2] = helikopter_4_B.travel;
    helikopter_4_B.TmpSignalConversionAtToWorkspac[3] = helikopter_4_B.Gain_d;
    helikopter_4_B.TmpSignalConversionAtToWorkspac[4] = helikopter_4_B.Gain_i;
    helikopter_4_B.TmpSignalConversionAtToWorkspac[5] = helikopter_4_B.Gain_b;
    helikopter_4_B.TmpSignalConversionAtToWorkspac[6] = helikopter_4_B.elevation;
    helikopter_4_B.TmpSignalConversionAtToWorkspac[7] = helikopter_4_B.Gain_dg;
  }

  /* Gain: '<S2>/Gain1' */
  rtb_Gain1_idx_4 = helikopter_4_P.Gain1_Gain_f * helikopter_4_B.elevation;
  rtb_Gain1_idx_5 = helikopter_4_P.Gain1_Gain_f * helikopter_4_B.Gain_dg;

  /* Sum: '<Root>/Sum1' incorporates:
   *  Constant: '<Root>/Vd_bias'
   *  Gain: '<S2>/Gain1'
   *  Gain: '<S6>/K_pd'
   *  Gain: '<S6>/K_pp'
   *  Sum: '<S6>/Sum2'
   *  Sum: '<S6>/Sum3'
   */
  helikopter_4_B.Sum1 = ((helikopter_4_B.u[0] - helikopter_4_P.Gain1_Gain_f *
    helikopter_4_B.Gain_i) * helikopter_4_P.K_pp - helikopter_4_P.Gain1_Gain_f *
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
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5181835955883618,
      0.51560932592839881, 0.50946586620907863, 0.50194000280169393,
      0.49916939222800688, -0.088716229908233926, 0.48558062273695429,
      0.03917630845735301, -0.011950840974435212, 0.516730423470277,
      -0.42755564960367404, -0.29632304342375176, -0.11053897235773937,
      -0.029583194858969552, -0.034068500551676269, -0.0759990330430028,
      -0.10710769475707842, -0.11367380931340236, -0.0987574668485768,
      -0.059394482552091635, -0.026429999470761616, -0.00925654409174438,
      -0.010294910641441777, -0.027911224107434549, -0.058631317801534126,
      -0.088749598397911175, -0.10889156170889018, -0.10929262154304062,
      -0.09645084847708342, -0.076521278735038376, -0.055384543652362672,
      -0.038401801036652021, -0.026782721132520733, -0.020563193201725814,
      -0.019953922012083364, -0.024732416846026471, -0.033609027494220338,
      -0.047431143270318145, -0.061111784710323372, -0.074274206750614385,
      -0.083589381422323042, -0.0882452857078776, -0.088238282925076572,
      -0.084573906361825463, -0.078988478824774028, -0.072539310655730091,
      -0.066157475105755981, -0.060237454130291489, -0.054863647553520051,
      -0.050071097903478634, -0.045727634472981581, -0.041713414565977977,
      -0.037978382147588877, -0.034472294638731804, -0.031178846844594158,
      -0.028043854871792055, -0.025119024603893605, -0.0223564279115352,
      -0.019831354633399623, -0.017412035340136318, -0.015235845257351782,
      -0.01325678013829352, -0.011435128387166722, -0.0097910512835827756,
      -0.00833911230670386, -0.00701472957488352, -0.005862455851704495,
      -0.0048506789860391893, -0.003981493390691278, -0.0032164434948575597,
      -0.0025531528718761329, -0.0019860946825908414, -0.0014937608853218098,
      -0.0010720725621434679, -0.00072030920070218342, -0.00045883292499082005,
      -0.00021662854903731259, -3.3165494145183395E-5, 8.83290297012384E-5,
      0.00019874784797924608, 0.00029728073945495172, 0.00031976836699097912,
      0.00035395949955027261, 0.00036795564329340514, 0.00039516196427605093,
      0.00038192384987935249, 0.00036629080299108234, 0.00036228717459441654,
      0.00032465684627129842, 0.00034757309241258589, 0.00030372930224737191,
      0.000283485336104314, 0.00025091960617664919, 0.0002432898024216129,
      0.00020406147184498191, 0.00015477089117772085, 0.00016122919481121104,
      7.119860935113817E-5, 0.00016091207369171046, 0.00019791215279090798,
      0.00019791215279090798, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      -1.9922682747042949E-15, 1.4610728030779248E-15, -0.14252494842275942,
      -0.33388750959003832, -0.53310859900500229, -0.70998092928305667,
      -0.81517897538186457, -0.642420324663873, -0.38618019371653539,
      0.11516615102543248, 0.70610225196191745, 1.2834477427951574,
      1.6149028520186779, 1.4694724851768011, 1.0183871925416808,
      0.47072639390109661, -0.071113710857005075, -0.54752311956066646,
      -0.89546101644131459, -1.0747745692983293, -1.0930638416456893,
      -0.99728752045750912, -0.83629841839455143, -0.63648748905515384,
      -0.40644684787545071, -0.15380316002589714, 0.10207005482461028,
      0.3284681832795025, 0.49366083813101092, 0.583700101699854,
      0.60801703928425477, 0.58912201090460359, 0.54913184854662966,
      0.50156805321961506, 0.45104681081976056, 0.39684005854005339,
      0.33651394947056257, 0.2686558314975056, 0.19454090224331341,
      0.11802745473914926, 0.045733706347877941, -0.015557589135845703,
      -0.060725908280623388, -0.088024893513631727, -0.099207244642779274,
      -0.0983708999212737, -0.090235897201733278, -0.078828204193638271,
      -0.066884533465129445, -0.05589174206695123, -0.046437779459302712,
      -0.038571408104629459, -0.032104449615341264, -0.026789724275221735,
      -0.022392792521255381, -0.018721466544009729, -0.015625086800776754,
      -0.013000027598884059, -0.010759247611640597, -0.0088443997621022123,
      -0.0071945835058289795, -0.0057893406400111348, -0.0045890364615269411,
      -0.0035614567649600755, -0.0026893663023910917, -0.0019554601966249867,
      -0.0013374234588011701, -0.00082754238967922228, -0.00040932558713426649,
      -6.9550627875046778E-5, 0.00020658733401035905, 0.00042622411040777971,
      0.00059548131033986147, 0.00072110930150267336, 0.00080652370382576955,
      0.00085501405446925091, 0.00087170679375433492, 0.00086984150421499592,
      0.00084759526431605, 0.00081001574623469847, 0.0007670934216628297,
      0.000717766127465034, 0.00065830593889581611, 0.00060292367046920414,
      0.00054904936494316544, 0.00049742963830371107, 0.00044171669318016944,
      0.00038849154260035942, 0.00033911202187329709, 0.00028944503714457208,
      0.00024624891352342012, 0.00019463457126314739, 0.00014719446530950671,
      0.00010192235197896031, 6.1248222735831387E-5, 1.9138256047491883E-5,
      -1.7804355575131474E-5, -4.4208936789272579E-5, -7.3406526372398842E-5,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helikopter_4_DW.FromWorkspace_PWORK.TimePtr = (void *) pTimeValues0;
    helikopter_4_DW.FromWorkspace_PWORK.DataPtr = (void *) pDataValues0;
    helikopter_4_DW.FromWorkspace_IWORK.PrevIndex = 0;
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
      3.1415926535897931, 3.1378809300300832, 3.1263930263617,
      3.1038248420568233, 3.0678252720296326, 3.0168274368423984,
      2.9540299741431153, 2.87931170631968, 2.7947693222151551,
      2.7031120570364635, 2.60289890329178, 2.4984890836981033,
      2.3941736353193326, 2.2920232719080449, 2.192437857635428,
      2.0952629802032812, 2.0005081735097456, 1.9084534643738416,
      1.8194510310473002, 1.7336967164682562, 1.6510174588307569,
      1.570927524691607, 1.4928219458450025, 1.4161896293995497,
      1.3407775608422181, 1.2666820763256, 1.1943014182674558,
      1.1241887875317884, 1.056835109509455, 0.99252107388139488,
      0.93126417983995036, 0.87284736647641081, 0.81690492673704485,
      0.763017016890854, 0.7107878850454169, 0.65990484266818961,
      0.61017787544384838, 0.5615561386312754, 0.51414018947979256,
      0.46815401525831091, 0.42390786517515822, 0.38174198134778181,
      0.34197012973281454, 0.30483447948143649, 0.27048025470931808,
      0.23895393763326028, 0.21021653332296661, 0.18416537371097361,
      0.16065691507097074, 0.13952582408272973, 0.12059932727423635,
      0.10370627514113236, 0.088681913544530214, 0.075370031803567092,
      0.063623522364257445, 0.053304229843926035, 0.044282227563698366,
      0.03643543776499502, 0.0296490893897684, 0.023815775820677853,
      0.018834738558831648, 0.014612045814871607, 0.011060676120496687,
      0.00810027182513219, 0.0056569870640322435, 0.0036634994079337417,
      0.0020586060704889135, 0.00078709534167088757, -0.00020049912450925662,
      -0.00094811469994321722, -0.0014945989387532945, -0.0018741442135229995,
      -0.0021166725588380756, -0.0022482815827645258, -0.0022916477713913778,
      -0.0022663403966001073, -0.0021888989620960785, -0.0020733321152204484,
      -0.0019313603049989619, -0.0017725086688264394, -0.0016045003465609102,
      -0.0014336552242303253, -0.0012647500044262548, -0.0011013722038735914,
      -0.00094612339212368714, -0.00080097494941025048, -0.000667213595240013,
      -0.0005455463998052, -0.00043633322987356209, -0.00033948650596907741,
      -0.00025496517891270123, -0.00018243653649059889, -0.00012140200611232992,
      -7.1180536654424464E-5, -3.1119138181364569E-5, -4.2454913013653503E-7,
      2.1931354763417137E-5, 3.6793288056389138E-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, -8.4785878377133276E-23, -1.7692965031688646E-22,
      -0.014846894238838998, -0.045951614673531788, -0.090272737219507418,
      -0.14399828010876345, -0.2039913407489376, -0.2511898507971323,
      -0.29887307129374174, -0.33816953641810027, -0.36662906071476536,
      -0.40085261497873265, -0.41763927837470743, -0.41726179351508291,
      -0.4086014536451516, -0.3983416570904672, -0.38869950972858763,
      -0.37901922677414263, -0.36821883654361626, -0.356009733306166,
      -0.34301725831617547, -0.3307170305499984, -0.320359736556599,
      -0.31242231538641874, -0.30652926578181061, -0.3016482742293265,
      -0.29638193806647323, -0.28952263223257535, -0.2804505229426697,
      -0.26941471208933349, -0.25725614251224094, -0.24502757616577822,
      -0.23366725345415787, -0.22376975895746409, -0.21555163938476349,
      -0.20891652738174854, -0.2035321695089089, -0.19890786889736478,
      -0.19448694725029222, -0.18966379660593141, -0.18394469688592668,
      -0.1769846003326109, -0.16866353530950562, -0.15908740645986916,
      -0.1485426010055122, -0.13741689908847363, -0.12610526830423111,
      -0.11494961724117468, -0.104204638447972, -0.094033834560011423,
      -0.084524363952964157, -0.075705987233973482, -0.067572208532415912,
      -0.060097446386408633, -0.05324752696385248, -0.046986037757238579,
      -0.0412771700813256, -0.036088009120910688, -0.031387159194813392,
      -0.027145393500906462, -0.023333254276362195, -0.01992414904738481,
      -0.016890770975840179, -0.014205478777499675, -0.011841617181457998,
      -0.0097731390443997826, -0.007973950624394004, -0.0064195733497793128,
      -0.0050860429152721036, -0.003950377864720577, -0.0029904623017358422,
      -0.002185936955240309, -0.0015181810990788202, -0.00097011338126030532,
      -0.00052643609570580062, -0.00017346475450740687, 0.00010122949916508239,
      0.0003097657380161149, 0.00046226738750251996, 0.00056788724088594682,
      0.00063540654469008983, 0.00067203328906211739, 0.00068338048932233854,
      0.00067562087921628261, 0.00065351120221065269, 0.0006209952469996171,
      0.0005805937708537463, 0.00053504541668095, 0.000486668781739252,
      0.00043685267972655157, 0.00038738689561793868, 0.00033808530822550476,
      0.00029011456968840936, 0.00024413812151307591, 0.00020088587783162183,
      0.00016024559389223961, 0.00012277835620491212, 8.942361557421469E-5,
      5.944773317188801E-5, 3.4524103049327084E-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 5.414806833785657E-23, 0.10493217810664326,
      0.21983628441780836, 0.313244766889668, 0.37971161810959081,
      0.42400803986259672, 0.33358104281428785, 0.33700679114151433,
      0.27773240726695775, 0.20114104837602229, 0.24187901076810089,
      0.11864172566642024, -0.0026679188175979984, -0.061207974615892141,
      -0.072512323594099315, -0.068147014993916288, -0.068416542797300972,
      -0.076333033229679756, -0.086289278742428172, -0.09182585110152966,
      -0.086933312108896654, -0.073201398254402819, -0.056098661345480104,
      -0.041649823912963221, -0.034496984129007313, -0.037220452663047326,
      -0.048478953905904357, -0.064118203611973235, -0.077996896279100211,
      -0.085932126130997738, -0.086426836563514964, -0.080290422154843083,
      -0.069951711019786145, -0.058082530449260192, -0.046894437534050371,
      -0.038054585034001082, -0.032682790594671315, -0.031245385748932736,
      -0.034088186681177371, -0.040420412564083436, -0.049191304216435042,
      -0.058810109575463718, -0.067680421363383717, -0.074526657645472311,
      -0.078632211985865633, -0.079946286208713566, -0.0788438811116459,
      -0.075941406353590243, -0.0718833574140604, -0.067209306363093191,
      -0.062324918707103667, -0.057486441383757116, -0.0528287640628445,
      -0.048412614335429177, -0.044253814304304827, -0.0403480963842983,
      -0.036675007807213812, -0.033223811913973417, -0.029979179895774456,
      -0.026942744094638368, -0.024094253742008034, -0.021438757692194468,
      -0.018978619682460923, -0.016706870947219451, -0.01461921347292248,
      -0.012715976600788374, -0.010985744924222074, -0.00942488380487596,
      -0.0080264468404593146, -0.0067843165851452511, -0.0056860778821294896,
      -0.0047194433598909177, -0.0038735333097604835, -0.0031357416036471961,
      -0.0024946666316434216, -0.0019414340728466789, -0.0014738544913698527,
      -0.0010778234146507732, -0.00074648078504087, -0.00047720065209859718,
      -0.00025886384061779636, -8.0197677680122608E-5, 5.4841960654421814E-5,
      0.00015626275287703128, 0.00022981035283368634, 0.0002855422031376533,
      0.00032191837131979742, 0.00034190757960833124, 0.000352081183100295,
      0.00034960526994874219, 0.00034844479026932584, 0.00033903885884185674,
      0.00032494397623069492, 0.00030569033930487425, 0.00028723000541595309,
      0.00026480412636643277, 0.00023573856782823873, 0.0002118580883344869,
      0.00017615069879349887, 0.00016226339187423787, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.41972871242657306, 0.45961642524466034,
      0.37363392988743865, 0.26586740487969113, 0.17718568701202364,
      -0.36170798819323563, 0.013702993308906048, -0.23709753549822662,
      -0.30636543556374185, 0.16295184956831443, -0.49294914040672261,
      -0.485238577936073, -0.23416022319317659, -0.045217395912828619,
      0.017461234400732006, -0.0010781112135386333, -0.03166596172951517,
      -0.039824982050993642, -0.022146289436405956, 0.019570155970531997,
      0.05492765541797532, 0.068410947635690875, 0.057795349730067544,
      0.02861135913582365, -0.010893874136160069, -0.045034004971428142,
      -0.062556998824275525, -0.055514770668507935, -0.031740919407590051,
      -0.0019788417300689097, 0.024545657634687506, 0.041354844540227732,
      0.047476722282103881, 0.044752371660839264, 0.035359410000197156,
      0.021487177757319079, 0.005749619382954309, -0.011371203728978504,
      -0.025328903531624262, -0.035083566609406418, -0.038475221436114725,
      -0.03548124715167994, -0.027384945128354408, -0.016422217361573264,
      -0.0052562968913917275, 0.004409620388270607, 0.011609899032222698,
      0.016232195758119346, 0.018696204203868787, 0.019537550623958119,
      0.019353909293386214, 0.018630709283650432, 0.017664598909661315,
      0.016635200124497415, 0.015622871680026108, 0.014692354308337961,
      0.013804783572961601, 0.012978528072795838, 0.012145743204544354,
      0.011393961410521326, 0.010621984199254267, 0.0098405520389341829,
      0.009086994940965893, 0.0083506298971878917, 0.00761294748853642,
      0.0069209267062652, 0.0062434444773844626, 0.0055937478576665841,
      0.00496852102125625, 0.0043929548120630463, 0.00386653808895429,
      0.0033836402005217345, 0.002951166824453151, 0.0025642998880150984,
      0.0022129302351869703, 0.0018703183259073046, 0.0015841243068763174,
      0.0013253705184396139, 0.0010771205317690912, 0.00087334724592320325,
      0.000714664651750695, 0.00054015855333817769, 0.00040568316889043793,
      0.00029419039982662027, 0.00022292740121586794, 0.0001455046727285764,
      7.9956833154135179E-5, 4.0694413967855065E-5, -9.9036526062111568E-6,
      -4.6419187176654185E-6, -3.7623725709876411E-5, -5.6379530444647194E-5,
      -7.70145477032828E-5, -7.3841335555684747E-5, -8.9703516198081314E-5,
      -0.00011626223415277624, -9.5521917975007189E-5, -0.00014282955816395218,
      -5.5549227677044094E-5, 1.2071854970196994E-5, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -3.1135087614589863E-17,
      -6.226427792322419E-17, -6.2289461451867394E-17, -0.0022269523191056776,
      -0.0091141588957792362, -0.022574589557751542, -0.043620955841531979,
      -0.071790174082869085, -0.10227332790172097, -0.13004803732263592,
      -0.14748157753008539, -0.14749188441560135, -0.12514134394965637,
      -0.080841020843412034, -0.022699937433627883, 0.038081315957852309,
      0.09137704242856752, 0.12964266498762178, 0.1483590668758307,
      0.14637912326966024, 0.12578270249981116, 0.090969090595622,
      0.04731090943377201, 7.8718734413047453E-5, -0.046029774266495641,
      -0.086963105995456272, -0.11934706394466739, -0.1406813892687622,
      -0.14968502002395564, -0.14652614578722928, -0.13269784758274997,
      -0.11053688666264858, -0.082637730683707625, -0.051406039711856921,
      -0.018854056109479485, 0.013410757381860114, 0.044104588041763766,
      0.072173448413077476, 0.096733706870558725, 0.11706589217976709,
      0.13264774597212012, 0.14321957091276158, 0.14883273125718074,
      0.14985595340309846, 0.1469224696254928, 0.14083074432282175,
      0.13242924144665022, 0.12251769801570038, 0.1117831428543585,
      0.10077281661646413, 0.089895151861147882, 0.079436737730976842,
      0.0695856371338833, 0.060454480633776911, 0.052100248236679858,
      0.044540085295809682, 0.037763373796708778, 0.031740759358373942,
      0.02643062038281668, 0.021783953542242253, 0.017747781222047072,
      0.014267862340524576, 0.011290155650288045, 0.0087622365888285914,
      0.0066342408487457486, 0.0048593127485083085, 0.0033939025944963118,
      0.0021980209757481027, 0.0012351496838092033, 0.00047225642480999192,
      -0.00012029946180948312, -0.00056886745631783381, -0.00089675402138320189,
      -0.0011244759957032015, -0.0012699883620231093, -0.0013489507664579002,
      -0.0013749694070263002, -0.0013597856130742841, -0.0013133225971221265,
      -0.0012439850089487844, -0.0011588046562038276, -0.0010634962911668033,
      -0.00096269359889420932, -0.00086018841984503527, -0.00075884676572435142,
      -0.000660821184746523, -0.00056767268020021307, -0.00048058414744787617,
      -0.00040032768190236456, -0.00032732758009758746, -0.0002617998050088963,
      -0.00020369184097954966, -0.000152979080828289, -0.00010946191217907698,
      -7.2841200804554709E-5, -4.2708321415617556E-5, -1.8671482860601459E-5,
      -2.5472947808179857E-7, 1.3158812841172463E-5, 2.2075972754139457E-5, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.2451676711132704E-16,
      -1.2458257772324741E-19, -0.00890780927642246, -0.027548826306694238,
      -0.053841722647889216, -0.084185465135121748, -0.11267687296534837,
      -0.12193261527540759, -0.11109883768365984, -0.069734160829797789,
      -4.1227542063770149E-5, 0.089402161863779847, 0.17720129242497734,
      0.23256433363913659, 0.24312501356592078, 0.21318290588286085,
      0.15306249023621704, 0.074865607552835667, -0.00791977442468177,
      -0.082385683079396324, -0.13925444761675659, -0.17463272464739998,
      -0.18892876279743584, -0.18443397200363473, -0.16373332691584255,
      -0.1295358317968445, -0.085337301296379217, -0.036014523020773785,
      0.012635496946905482, 0.05531319281791721, 0.088643843680405648,
      0.11159662391576385, 0.12492676388740277, 0.13020793440950976,
      0.12905925396535839, 0.1227753226396146, 0.11227544148525485,
      0.098241033829925, 0.0813287412368335, 0.06232741516941203,
      0.042287299762565948, 0.022452641377676597, 0.0040928885836708836,
      -0.011733935110422616, -0.024366901210684317, -0.033606011504686138,
      -0.039646173723799305, -0.042938220645367509, -0.044041304951577487,
      -0.043510659021264976, -0.041833656520684159, -0.039404402388374206,
      -0.036524626000425525, -0.033416929588388218, -0.030240651763480678,
      -0.02710684599640361, -0.024090457753339359, -0.02124055590222907,
      -0.018586667362297714, -0.016144689280780715, -0.013919675526089987,
      -0.011910826760946128, -0.010111676245837817, -0.00851198296033137,
      -0.007099712400949758, -0.0058616406160479891, -0.0047835264749928347,
      -0.0038514851677555972, -0.0030515730359968458, -0.0023702235464779,
      -0.0017942719780334024, -0.0013115462602614728, -0.00091088789727999866,
      -0.00058204946527963173, -0.00031584961773916314, -0.00010407456227359984,
      6.0735175808064809E-5, 0.00018585206380862963, 0.0002773503526933681,
      0.00034072141097982764, 0.00038123346014809668, 0.00040321076909037638,
      0.000410020716196696, 0.0004053666164827353, 0.00039210232391131408,
      0.00037259401818523941, 0.00034835413100934784, 0.00032102586218204669,
      0.00029200040721910831, 0.00026211110035476474, 0.00023243185611738655,
      0.00020285104060504264, 0.000174068674596848, 0.0001464828454980891,
      0.00012053151755574861, 9.6147354220064372E-5, 7.3667013530078653E-5,
      5.3654169277017044E-5, 3.566863965186797E-5, 2.5929053936327696E-5, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0 } ;

    helikopter_4_DW.FromWorkspace1_PWORK.TimePtr = (void *) pTimeValues0;
    helikopter_4_DW.FromWorkspace1_PWORK.DataPtr = (void *) pDataValues0;
    helikopter_4_DW.FromWorkspace1_IWORK.PrevIndex = 0;
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
  helikopter_4_M->Sizes.checksums[0] = (33496848U);
  helikopter_4_M->Sizes.checksums[1] = (2747364283U);
  helikopter_4_M->Sizes.checksums[2] = (2217839343U);
  helikopter_4_M->Sizes.checksums[3] = (3845084763U);

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

    helikopter_4_B.TravelCounttorad = 0.0;
    helikopter_4_B.Gain = 0.0;
    helikopter_4_B.travel = 0.0;
    helikopter_4_B.Gain_d = 0.0;
    helikopter_4_B.PitchCounttorad = 0.0;
    helikopter_4_B.Gain_i = 0.0;
    helikopter_4_B.Gain_b = 0.0;
    helikopter_4_B.ElevationCounttorad = 0.0;
    helikopter_4_B.Gain_e = 0.0;
    helikopter_4_B.elevation = 0.0;
    helikopter_4_B.Gain_dg = 0.0;
    helikopter_4_B.u[0] = 0.0;
    helikopter_4_B.u[1] = 0.0;
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
  helikopter_4_M->Sizes.numBlocks = (61);/* Number of blocks */
  helikopter_4_M->Sizes.numBlockIO = (19);/* Number of block outputs */
  helikopter_4_M->Sizes.numBlockPrms = (154);/* Sum of parameter "widths" */
  return helikopter_4_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
