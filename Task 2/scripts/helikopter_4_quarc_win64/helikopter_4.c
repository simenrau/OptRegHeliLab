/*
 * helikopter_4.c
 *
 * Code generation for model "helikopter_4".
 *
 * Model version              : 1.180
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Tue Feb 27 18:42:11 2018
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
              pDataValues += 81;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 6; ++elIdx) {
              (&rtb_Sum3_j[0])[elIdx] = pDataValues[currTimeIndex + 1];
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
          for (elIdx = 0; elIdx < 6; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&rtb_Sum3_j[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1, f2);
            pDataValues += 81;
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

  /* Gain: '<S3>/Gain1' */
  tmp[0] = helikopter_4_B.travel;
  tmp[1] = helikopter_4_B.Gain_d;
  tmp[2] = helikopter_4_B.Gain_i;
  tmp[3] = helikopter_4_B.Gain_b;
  tmp[4] = helikopter_4_B.Sum;
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
    helikopter_4_B.TmpSignalConversionAtToWorkspac[6] = helikopter_4_B.Sum;
    helikopter_4_B.TmpSignalConversionAtToWorkspac[7] = helikopter_4_B.Gain_dg;
  }

  /* Gain: '<S2>/Gain1' */
  rtb_Gain1_idx_4 = helikopter_4_P.Gain1_Gain_f * helikopter_4_B.Sum;
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
      19.75, 20.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.52359877559829882,
      0.52359877559829882, 0.52359877559829882, 0.518830512689434,
      0.45655353970611007, 0.44135188285902188, 0.44604711319292312,
      0.43583034891715616, 0.37046117404703838, 0.2163450194773911,
      -0.018747175032251127, -0.19982634788849066, -0.33346563573471683,
      -0.42623504844971843, -0.48439856205140897, -0.51382434064353222,
      -0.51999645090359758, -0.507780639115384, -0.48163924748243364,
      -0.44538596763823918, -0.40237516913386423, -0.3554401169420448,
      -0.30694793522510139, -0.25884341644544356, -0.21268510500933596,
      -0.16969568208035571, -0.13080376051993764, -0.096680553783877368,
      -0.0677640869413091, -0.04428646767589977, -0.026261666682917868,
      -0.01347598812748753, -0.0054432637182297659, -0.0013516031239396659,
      -3.0579260831148543E-5, -3.3452296738592063E-5, -8.6923521043381544E-5,
      -8.6923521043381544E-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.285742826072486, 0.3022528665842143, 0.31764699488765646,
      0.33126944844502787, 0.34229690536613228, 0.34970651605462505,
      0.352240718820672, 0.34835239039617172, 0.33615632280092717,
      0.31335722882270833, 0.2771866464842348, 0.22429533078898875,
      0.15066502261241613, 0.051472439069692705, -0.079048986205959029,
      0.00019557303580803698, 9.1137987698735061E-5, 3.9939848758152125E-5,
      2.5666764273285807E-5, 3.6570772412832165E-5, 6.0276655988683706E-5,
      8.914636322980202E-5, 0.00011692258859541957, 0.00014205892957190662,
      0.00015985171837632763, 0.00017312903258281476, 0.0001759126049504339,
      0.00017436285113004763, 0.00016652450272873622, 0.00015559559740187012,
      0.00013985611228204846, 0.00011890073225372658, 9.93764738934852E-5,
      7.8065619279121167E-5, 5.7092284918891867E-5, 3.8476085387853542E-5,
      2.0634468627296696E-5, 9.58305098213413E-6, 2.9438402929562177E-8, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

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
      19.75, 20.0 } ;

    static real_T pDataValues0[] = { 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1415926535897931, 3.1415926535897931,
      3.1415926535897931, 3.1378421413625261, 3.1262155534579983,
      3.1033093000299643, 3.0666274151911783, 3.0144539223941584,
      2.9456562771175667, 2.8595419181257866, 2.7561377096658859,
      2.6357915909022895, 2.4988445367776211, 2.345639550310564,
      2.1769190260798039, 1.9944997661292598, 1.80183963215139,
      1.6034897505814103, 1.4043676089693413, 1.2091702306052616,
      1.0219934235060744, 0.84613250234711057, 0.6840175418677602,
      0.53723729033420653, 0.40661702373053638, 0.2923242688317208,
      0.19398557487487725, 0.1108032757054615, 0.041665560301247066,
      -0.014753839747948985, -0.059908126569292917, -0.09530255742655841,
      -0.12243190258401694, -0.14272949185837377, -0.15752694727008579,
      -0.16802363980195029, -0.1752651726351275, -0.1801304827320265,
      -0.18332764379809646, -0.18539896684140728, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.015002048909068423,
      -0.046506351618112111, -0.091625013712135384, -0.14672753935514371,
      -0.20869397118807925, -0.27519058110636674, -0.34445743596712158,
      -0.41361683383960285, -0.48138447505438448, -0.547788216498674,
      -0.61281994586822852, -0.67488209692303991, -0.72967703980217669,
      -0.77064053591147919, -0.79339952627991939, -0.7964885664482757,
      -0.78078951345631853, -0.74870722839674975, -0.70344368463585527,
      -0.64845984191740147, -0.5871210061342147, -0.52248106641468051,
      -0.45717101959526235, -0.39335477582737416, -0.33272919667766304,
      -0.27655086161685771, -0.22567760019678421, -0.18061714728537573,
      -0.14157772342906194, -0.10851738062983415, -0.08119035709742739,
      -0.05918982164684803, -0.041986770127458009, -0.028966131332708804,
      -0.019461240387595974, -0.012788644264279794, -0.00828529217324328,
      -0.0053491890494525678, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.10602875205865551, 0.22266037932317656, 0.31888147181640641,
      0.38944360631144165, 0.43795507377677839, 0.46997264230390062,
      0.48955167553642559, 0.48879221058362549, 0.47895580607117128,
      0.46931628339830855, 0.45961942605733425, 0.43863158067960462,
      0.3872697288093932, 0.2895143455907212, 0.16085185173753255,
      0.021832156133812473, -0.11095491071447786, -0.22674533783181425,
      -0.31990544010418137, -0.38860480072845421, -0.43351946459065,
      -0.45685040644455038, -0.46158646749707721, -0.45102883804691596,
      -0.42847843911450578, -0.39704701639984735, -0.35955278203789637,
      -0.31847007154512275, -0.27591573775444334, -0.23365787639277547,
      -0.19313696550256904, -0.15549138204962645, -0.12158459789443167,
      -0.092024902117153723, -0.067176939061619728, -0.047159360980284226,
      -0.031827972644770469, -0.020751255513029844, -0.013205193429379328, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.42411500823462206,
      0.46652650905808424, 0.38488436997291947, 0.282248537980141,
      0.19404586986134689, 0.12807027410848892, 0.078316132930099966,
      -0.0030378598112003651, -0.039345618049817038, -0.038558090691450671,
      -0.038787429363897286, -0.083951381510918577, -0.20544740748084578,
      -0.39102153287468777, -0.5146499754127547, -0.55607878241488029,
      -0.53114826739316134, -0.46316170846934551, -0.3726404090894686,
      -0.27479744249709137, -0.17965865544878318, -0.093323767415601624,
      -0.018944244210107434, 0.042230517800645087, 0.090201595729640871,
      0.12572569085863389, 0.14997693744780377, 0.1643308419710946,
      0.17021733516271745, 0.16903144544667165, 0.16208364356082564,
      0.1505823338117703, 0.13562713662077913, 0.11823878310911186,
      0.099391852222135937, 0.080070312325342, 0.061325553342055049,
      0.044306868526962494, 0.030184248334602069, 0.01975653374696924, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0044647316573825941,
      0.012535981440797888, 0.023482891641332385, 0.036673284713674356,
      0.051547548482380828, 0.067594390548613034, 0.084327852884822987,
      0.10126479338459858, 0.11790231860186939, 0.13369440681854297,
      0.1480272906042101, 0.16019259288049867, 0.1693577841503433,
      0.17453292519943295, 0.17453292519943295, 0.17180890407187632,
      0.16704023530102549, 0.16077984365690118, 0.15347494729042116,
      0.14548466137674104, 0.13709484271281791, 0.12853067379278946,
      0.11996726710082717, 0.11153863997461766, 0.10334517876460973,
      0.0954599967486345, 0.0879340904529084, 0.080800822701465311,
      0.074079503669911437, 0.067778432722745044, 0.061897322519282282,
      0.056429309679333794, 0.051362707142413228, 0.046682267051284448,
      0.042370286750789507, 0.038407492291576174, 0.034773683130257388,
      0.031448358927384019, 0.028411027436293763, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.017858926629530376, 0.032284999133661174,
      0.043787640802138, 0.05276157228936787, 0.059497055074825896,
      0.064187368264928837, 0.066933849344839841, 0.0677477619991023,
      0.066550100869083237, 0.063168352866694283, 0.057331535142668547,
      0.048661209105154278, 0.036660765079378585, 0.020700564196358566,
      2.914479494641829E-17, -0.010896084510226535, -0.019074675083403287,
      -0.025041566576497351, -0.029219585465920028, -0.031961143654720539,
      -0.033559274655692431, -0.034256675680113771, -0.034253626767849238,
      -0.03371450850483803, -0.0327738448400317, -0.031540728063900951,
      -0.030103625182904419, -0.028533071005772341, -0.026885276126215489,
      -0.025204283788665582, -0.023524440813851022, -0.021872051359793972,
      -0.020266410147682278, -0.018721760364515127, -0.017247921201979741,
      -0.015851177836853345, -0.014535236645275146, -0.013301296811493488,
      -0.012149325964361019, -0.011077516906232265, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

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

  rtmSetTFinal(helikopter_4_M, 20.0);
  helikopter_4_M->Timing.stepSize0 = 0.002;
  helikopter_4_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  helikopter_4_M->Sizes.checksums[0] = (116918794U);
  helikopter_4_M->Sizes.checksums[1] = (3353328249U);
  helikopter_4_M->Sizes.checksums[2] = (1457290355U);
  helikopter_4_M->Sizes.checksums[3] = (1694034377U);

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
    helikopter_4_B.Sum = 0.0;
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
