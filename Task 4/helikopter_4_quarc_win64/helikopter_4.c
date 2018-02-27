/*
 * helikopter_4.c
 *
 * Code generation for model "helikopter_4".
 *
 * Model version              : 1.180
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Tue Feb 27 21:35:31 2018
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
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52248860413619336,
      0.521497892572679, 0.52134445347608127, 0.51883418312776952,
      0.52083820204572451, 0.52290999110105507, 0.41580262525738171,
      0.52086875543565869, -0.41173345982497761, -0.34511451978488134,
      0.50949453527482824, 0.509581753113237, -0.511936436497421,
      -0.277266811138346, 0.51819302682394808, 0.52359142945560166,
      -0.5131261648910973, -0.504717910976412, -0.50996588025474376,
      -0.48930463874741659, -0.49086018632832351, -0.33047666924564384,
      -0.30211990294647378, -0.093458291176296948, -0.14452781656533487,
      -0.04230880085966001, -0.049929771538664125, -0.026276425828239816,
      0.0040396637634171895, 0.018072989196673238, 0.011917571520030432,
      0.012708432313354954, -0.012586768734338654, -0.020908628077584552,
      -0.047739267906344654, -0.069297498973638877, -0.090227012290735542,
      -0.10760344593762465, -0.1152108302381971, -0.13228748856374895,
      -0.10557842399683171, -0.076137348234214114, -0.09724418699120585,
      -0.0476174069925342, -0.019546206343223124, -0.032787073650868595,
      -0.018207818945682148, -0.01570028705822062, -0.0021935092364640638,
      -0.018669243768539814, -0.0005418730822463998, -0.012255421659634643,
      -0.0094763822771806484, -0.0072788047979677646, -0.0065476258626210411,
      -0.0070409943385991288, -0.0049033939326523163, -0.0045107383236791211,
      -0.0013694614840104402, -0.004359618166705808, -0.0034747134896490539,
      -0.0043251811151125159, -0.0042506386757921568, -0.0056806106269583247,
      -0.006915778557075339, -0.0077858027494898414, -0.0092451347018249145,
      -0.0084147582264944842, -0.0095612753125615488, -0.010590668452947526,
      -0.011564156323352354, -0.012296737190097884, -0.012749418264546562,
      -0.014088089861468839, -0.014389551733694097, -0.015258406978617182,
      -0.016612442707168971, -0.017546517114164158, -0.018880252640556285,
      -0.020054804928297313, -0.021619003208955542, -0.023129028522988017,
      -0.02504017671366977, -0.027100095455758647, -0.0298275335591384,
      -0.033195578618085191, -0.0372375824071986, -0.041430482453391225,
      -0.044956755939219308, -0.045161684541826869, -0.039579873079415906,
      -0.025910240014859871, -0.0062867509002845282, 0.013593617351267951,
      0.024850527441818419, 0.023656601168448761, 0.012900322538134532,
      0.0023651147468487662, -0.00010019240314433084, -0.000226668854361369,
      -0.000226668854361369, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      9.6866485232378125E-16, 3.6396429530241649E-14, -0.19161192118711431,
      -0.44946151573309789, -0.72015430998232988, -0.96245116857570523,
      -1.1065449884814327, -1.0235029485533131, -0.5171837609429647,
      0.48377050410684097, 1.805008522041919, 2.5619164926110467,
      2.7054298777784234, 2.2622335295681069, 0.69708869561065867,
      -1.1659409725963492, -2.5714410307761719, -3.5788825085626925,
      -3.4243286376377369, -1.8274988949621691, 0.23698592731954571,
      1.6640411393914731, 2.0695025794208135, 1.7550177144120989,
      1.212969905834534, 0.78935075597394422, 0.51270634640957946,
      0.3765342736605593, 0.30626210063231146, 0.26781701313334128,
      0.24294568780077302, 0.21475708389745635, 0.17175312336588797,
      0.11455404400901703, 0.044624150957088181, -0.027744867258064826,
      -0.093724814136322124, -0.14309008185793728, -0.16890026396502208,
      -0.16889203062622632, -0.14691108004521342, -0.10609574266844912,
      -0.062877681244812716, -0.029470187179496555, 0.0050468345673961726,
      0.023518141838557848, 0.025837334607626335, 0.027424952262394355,
      0.025265542250591036, 0.022219993948516172, 0.01533708818960094,
      0.01404583452221324, 0.0092583784836437772, 0.00805877698553544,
      0.0077534854487745358, 0.00735860835121665, 0.0069749678679099222,
      0.0070509357881973151, 0.0067551393004719167, 0.0063892364131557019,
      0.0050482681399383307, 0.0046303566451321989, 0.0043921747276469388,
      0.0046498726515566751, 0.0051587478450640222, 0.0063620856344028635,
      0.0083378840811767236, 0.010965385197766884, 0.014410826576340819,
      0.017898380671455567, 0.0218846324579868, 0.026516055605755463,
      0.031833164736719946, 0.037768263251943136, 0.044202550358139692,
      0.051429000568256694, 0.059182753016855452, 0.0675798493776633,
      0.076845264172019326, 0.086901793717157561, 0.0978541081441593,
      0.10966368876369742, 0.12242907258011967, 0.13612385421846937,
      0.15082131108010746, 0.16651382184328339, 0.18331183835855294,
      0.20130144534733696, 0.22048973465723862, 0.24057077688641992,
      0.26074180956139215, 0.27909652846677135, 0.29271597869371652,
      0.29817183840556427, 0.29348525666845227, 0.27975649404419056,
      0.26162494937404607, 0.24462305953954161, 0.23202637446375249, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0 } ;

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
      3.1415926535897927, 3.137850093466001, 3.1262552535003545,
      3.1034206655171181, 3.0668812856805006, 3.0149359202055361,
      2.9464440915661827, 2.861433090387155, 2.7595253472943404,
      2.6469774642725872, 2.5303087300562961, 2.4086762176056467,
      2.2791799975976015, 2.1462022176593081, 2.01375954183865,
      1.8792468613111044, 1.7385448470031897, 1.5953362393935517,
      1.4548579221449875, 1.3217940292371106, 1.1995649552847443,
      1.0905273261460278, 0.9950899416645822, 0.91278648422669073,
      0.84148570966206326, 0.77935011895524464, 0.72427097477848223,
      0.67453909816084445, 0.6287428575262376, 0.58562933212386181,
      0.54414680766107382, 0.503561497896569, 0.46337871664154623,
      0.42345746770490589, 0.38385190691693449, 0.3448564435896655,
      0.30693296637444406, 0.27065041104336796, 0.23661935216663491,
      0.20539381001926466, 0.17752082735836205, 0.15324324437024803,
      0.13245181985098595, 0.11510411656577914, 0.10083102484096211,
      0.0890382886720941, 0.07927872581060412, 0.071135767741403,
      0.064268500400743792, 0.058331181910380475, 0.053164637583073986,
      0.048568516718100738, 0.044454487458784296, 0.0407660847276008,
      0.037452176186729087, 0.034468359551834644, 0.031784793816421186,
      0.029366774037984676, 0.027182339111458349, 0.02518432925517218,
      0.02335195528853852, 0.021670276418271447, 0.020135238356575359,
      0.018746367033024593, 0.017514610210827918, 0.01646100128283453,
      0.015611597223224362, 0.014999151225118846, 0.014645861620489326,
      0.014574451178973508, 0.014810573444088815, 0.015382546721532355,
      0.016319008391304689, 0.017646135816161736, 0.019393650126365385,
      0.021588485269505124, 0.02425750151332836, 0.027431796524409715,
      0.031144143233934881, 0.0354307383039199, 0.040329897511055257,
      0.045884022994643087, 0.052138982822378992, 0.0591460946818999,
      0.0669626880355677, 0.075655838196199754, 0.085306110986314618,
      0.096010791202511547, 0.10788255449614745, 0.121041173955952,
      0.13558578348847414, 0.151549216791137, 0.16883968801635732,
      0.18720573845281815, 0.20625612663886059, 0.22555592053458209,
      0.24475359747330255, 0.26367109615751627, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.6882561540465575E-25,
      -0.014970240495168053, -0.046379359862587642, -0.091338351932945067,
      -0.14615751934647112, -0.20778146189985822, -0.27396731455741585,
      -0.34004400471610968, -0.40763097237125917, -0.45019153208701279,
      -0.46667493686516387, -0.48653004980259773, -0.51798488003218157,
      -0.5319111197531734, -0.52977070328263343, -0.53805072211018168,
      -0.56280805723165883, -0.572834430438552, -0.56191326899425664,
      -0.53225557163150827, -0.4889162958094655, -0.43615051655486514,
      -0.3817495379257822, -0.3292138297515661, -0.28520309825850976,
      -0.24854236282727429, -0.22031657670704941, -0.1989275064705511,
      -0.183184962538428, -0.17245410160950292, -0.16593009785115143,
      -0.16234123905801984, -0.16073112502009113, -0.15968499574656136,
      -0.15842224315188538, -0.1559818533090761, -0.15169390886088582,
      -0.14513022132430445, -0.13612423550693226, -0.12490216858948103,
      -0.1114919306436104, -0.09711033195245615, -0.083165698077048308,
      -0.069390813140827151, -0.057092366899268138, -0.047170944675472073,
      -0.039038251445959953, -0.032571832276804484, -0.027469069362636867,
      -0.023749273961453266, -0.020666177309225972, -0.018384483459892988,
      -0.016456117037265759, -0.014753610924733995, -0.013255634163486825,
      -0.01193526653957777, -0.010734262941653843, -0.00967207911374604,
      -0.00873773970610531, -0.0079920394251446758, -0.0073294958665346575,
      -0.00672671548106829, -0.0061401522467843511, -0.0055554852942030749,
      -0.0049270272887867016, -0.0042144357119735555, -0.003397616238440667,
      -0.0024497839924220617, -0.0014131584185180792, -0.00028564176606327352,
      0.00094448906046122529, 0.0022878931097741655, 0.0037458466790893426,
      0.0053085096994281763, 0.0069900572408146053, 0.00877934057255896,
      0.010676064975292934, 0.012697180044325409, 0.014849386838100668,
      0.017146380279940088, 0.019596636828541417, 0.022216501934351322,
      0.025019839310943613, 0.028028447438083637, 0.031266373414671159,
      0.034772600642528254, 0.038601091160459469, 0.042818720864787704,
      0.047487053174543666, 0.052634477839218187, 0.0581784381300886,
      0.063853733210651534, 0.069161884900881293, 0.07346420174584338,
      0.076201552744169745, 0.077199175582886012, 0.076790707754881823,
      0.075669994736854956, 0.074493007133470457, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 1.0319798253948019E-29, 0.10580394233757916,
      0.22198765981730453, 0.31775298507082378, 0.38744093856654488,
      0.43553449035353436, 0.46777632867047164, 0.46700480982019449,
      0.47767887435511708, 0.300801781205462, 0.11649840957714681,
      0.14032835511336197, 0.22231072673304175, 0.098425343593826167,
      -0.015127646139051419, 0.058520010741760525, 0.17497538923775896,
      0.070862576521705092, -0.077186598043269938, -0.20960927799699833,
      -0.30630544923496084, -0.37292837529585049, -0.38448533995410272,
      -0.37130232076917913, -0.31105104147318113, -0.25910407644230288,
      -0.19948907621497702, -0.15116977944925758, -0.11126228806937136,
      -0.0758416266808161, -0.046109166895586536, -0.025364683281392159,
      -0.011379671080162069, -0.0073936421642607254, -0.0089246530646486435,
      -0.017247743367459434, -0.030305553694311031, -0.046389636683959595,
      -0.063650868162195479, -0.079313283004752155, -0.0947784401203224,
      -0.10164364688277555, -0.098555346452083942, -0.097355625780971913,
      -0.086920720973301877, -0.070120823056390752, -0.057478769681906891,
      -0.045702181012040016, -0.036064379413756863, -0.026290093219374415,
      -0.021790149631779717, -0.016126140695251825, -0.013628957387239659,
      -0.012032662976779029, -0.010587127636405408, -0.0093318540867510513,
      -0.0084882347389798451, -0.00750710962299773, -0.0066035540873011223,
      -0.0052703247855865349, -0.0046826048314945242, -0.0042602215486579834,
      -0.0041456049178069726, -0.0041322027229665414, -0.0044417011595856469,
      -0.00503632510965459, -0.0057729680766170823, -0.0066989163157278224,
      -0.0073264736449890717, -0.0079688570747740266, -0.0086940948663641084,
      -0.0094946667433602281, -0.010304259001615929, -0.011044305410481733,
      -0.011884535800488528, -0.012645971219939711, -0.013405323675455881,
      -0.014284469397171424, -0.015210975640681405, -0.016234272371817838,
      -0.0173174774756737, -0.018516205980579674, -0.019812917918151415,
      -0.021263693185359336, -0.022884424163445591, -0.024780674936944476,
      -0.027058308791358888, -0.029808595939526786, -0.03299398982090513,
      -0.036380031608980692, -0.0391825939687749, -0.040110818102517452,
      -0.037515988840223918, -0.030407132305542393, -0.019346551398297138,
      -0.0070508172087320473, 0.0028868946049905083, 0.0079207715851032454,
      0.0083184988618408163, 0.0067340263818944394, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.42321576935031663, 0.46473486991890162,
      0.38306130101407687, 0.27875181398288434, 0.19237420714795808,
      0.12896735326774905, -0.0030860754011087605, 0.042696258139690642,
      -0.70750837259862043, -0.7372134865132608, 0.0953197821448606,
      0.32792948647871911, -0.4955415325568624, -0.45421195893151034,
      0.29459062752324777, 0.4658215139839938, -0.41645125086421553,
      -0.59219669825990018, -0.52969071981491356, -0.3867846849518502,
      -0.26649170424355834, -0.046227858633008989, 0.052732076739694184,
      0.24100511718399209, 0.20778786012351305, 0.23846000090930344,
      0.19327718706287775, 0.15962996551954489, 0.14168264555422105,
      0.11892983914091822, 0.082977934456777508, 0.055940048804920359,
      0.015944115663605374, -0.0061240436015516722, -0.033292361211243154,
      -0.052231241307406404, -0.064336331958594284, -0.06904492591294345,
      -0.062649659370226735, -0.061860628462280978, -0.027460827049812547,
      0.012353201722766466, 0.0047988826844480963, 0.041739619230680106,
      0.06719959166764454, 0.050568213497935409, 0.047106354679467506,
      0.038551206393132635, 0.039097144777529791, 0.017999774350378785,
      0.022656035746111573, 0.0099887332320486664, 0.0063851776418425189,
      0.0057821413614944873, 0.0050210941986174188, 0.0033744773910848271,
      0.0039245004639284582, 0.0036142221427864314, 0.0053329172068583485,
      0.0023508798163680397, 0.001689533131346164, 0.0004584665234040434,
      5.3608779361724183E-5, -0.0012379937464764225, -0.002378495800275768,
      -0.0029465718678499736, -0.0037037929564429609, -0.0025102293170449924,
      -0.0025695337191398164, -0.0029009511663603306, -0.0032022875079844777,
      -0.0032383690330228063, -0.0029601856354632113, -0.0033609215600271785,
      -0.0030457416778047327, -0.0030374098220646837, -0.0035165828868621704,
      -0.0037060249740399238, -0.0040931869245457329, -0.0043328204154234606,
      -0.0047949140196238858, -0.0051868477502869867, -0.0058031010688316814,
      -0.0064829239123450257, -0.0075850030939955447, -0.0091105354176576319,
      -0.0110011485926716, -0.012741575525513353, -0.013544167152302281,
      -0.011210249439176843, -0.0037128965349702003, 0.010379317049174157,
      0.028435426138726096, 0.044242323628981027, 0.049182936758260368,
      0.039750847254890219, 0.020135507920450949, 0.0015909091069502834,
      -0.0063378899197855084, -0.007555374842102322, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.1138301989743378E-16,
      -9.7971743931788262E-17, -1.5693686213808799E-16, -0.0029939362685482379,
      -0.012262224653288968, -0.030419071781122468, -0.058883409375785345,
      -0.097046130020224844, -0.14074035080320335, -0.18007566687360502,
      -0.19981917181843686, -0.18360986007525568, -0.1283007565111591,
      -0.041677682934122946, 0.060641720468642711, 0.14892449768547916,
      0.19597122601896602, 0.18875056088686784, 0.12435297243495579,
      0.01960041861907481, -0.08946168217091599, -0.16786160918996423,
      -0.19926307286733808, -0.18785535517832505, -0.14876392961032406,
      -0.097557465730997439, -0.044494575859248264, 0.0048379636092601395,
      0.048415943984388428, 0.085809181406719837, 0.11728225117892076,
      0.14334231142047857, 0.16441040086287398, 0.1806553868813176,
      0.19206012081930879, 0.19848818321148609, 0.19987477706691009,
      0.19634889437491862, 0.18834565643522427, 0.17663620988139195,
      0.16227228610568253, 0.14644391686879726, 0.13037938949153749,
      0.11506034398806751, 0.10107341022498012, 0.088864248817966576,
      0.078495576694168581, 0.06973427706678352, 0.062365323839498249,
      0.056143784937531291, 0.050850359981509476, 0.0462426866278067,
      0.042211860902228115, 0.038610861793291637, 0.035376470525390953,
      0.032468530569082353, 0.029849796505379434, 0.027487409040396427,
      0.0253593862429532, 0.02343942743018446, 0.021703049729517317,
      0.020113404590106871, 0.01865440990610551, 0.017314519676503558,
      0.016090781109699765, 0.01498304324973068, 0.014000228987952352,
      0.01315928767960995, 0.0124811622641315, 0.011992123497784083,
      0.011709988460637958, 0.011652957635281532, 0.011841529315406675,
      0.012298273811460505, 0.013045937401259261, 0.01410518941465019,
      0.015499363786677722, 0.01724933149698283, 0.019375564867070805,
      0.021901426342684233, 0.024850919775176915, 0.028249800502694468,
      0.032124160563778266, 0.036501481735800881, 0.0414094678281724,
      0.046876704730955822, 0.052931887939529507, 0.059605074308890969,
      0.066928238420409056, 0.074934434321990567, 0.0836522459117081,
      0.093093844842111889, 0.10322886095483762, 0.11394921888081323,
      0.12503547134796256, 0.13615541128875638, 0.14691288222398033,
      0.15694144695798096, 0.16598959202903707, 0.17394890282460682, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.05415532702363E-17,
      1.923660117600047E-15, -0.011975745074193194, -0.0370731535389635,
      -0.07262738851133399, -0.11385735037865151, -0.15265088257775802,
      -0.174776883131914, -0.15734126428160675, -0.0789740197793273,
      0.064837246972724763, 0.22123641425638629, 0.3464922943081446,
      0.40927761361106263, 0.35313110886734583, 0.18818691333394735,
      -0.028882660528392657, -0.25759035380764816, -0.41901021526352394,
      -0.43624840315996322, -0.313599708076193, -0.1256058547094954,
      0.045630870756052072, 0.15636570227200389, 0.20482585551730659,
      0.2122515594869967, 0.1973301578740336, 0.17431192150051314,
      0.14957294968932558, 0.12589227908880374, 0.10424024096623113,
      0.084272357769581821, 0.064979944073774457, 0.045618935751964788,
      0.025712249568709256, 0.0055463754216960887, -0.014103530767965949,
      -0.032012951758777421, -0.046837786215329358, -0.057455695102837684,
      -0.0633134769475411, -0.064258109509039055, -0.061276182013879917,
      -0.055947735052349561, -0.048836645628054128, -0.041474688495191987,
      -0.03504519850954025, -0.029475812909141079, -0.024886155607867839,
      -0.021173699824087261, -0.018430693414811093, -0.016123302902314337,
      -0.014403996435745935, -0.012937565071602744, -0.011631759825234378,
      -0.010474936254811678, -0.0094495498599320334, -0.0085120911897729069,
      -0.0076798352510749623, -0.0069455108026685655, -0.0063585805576418077,
      -0.0058359787360054262, -0.0053595609184078144, -0.0048949542672151632,
      -0.0044309514398763435, -0.0039312570471133138, -0.0033637652333696083,
      -0.0027125016619137981, -0.0019561550653896692, -0.0011285401485844979,
      -0.00022812330142570368, 0.00075428672050056639, 0.0018269779842153256,
      0.0029906543591950233, 0.0042370080535637164, 0.0055766974881101272,
      0.0069998708412204241, 0.0085049334803519172, 0.010103445902453718,
      0.011797973729970709, 0.013595522910070223, 0.015497440244335199,
      0.017509284688090473, 0.019631944369486055, 0.021868947611133702,
      0.024220732834294711, 0.026692745477445854, 0.029292656446072356,
      0.032024783606326, 0.034871246358870171, 0.037766395721615234,
      0.040540064450902878, 0.042881431703902441, 0.044345009868597247,
      0.044479759763175371, 0.04302988374089578, 0.040114258936002449,
      0.036192580284224417, 0.031837243182279033, 0.01350358288489446, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0 } ;

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
  helikopter_4_M->Sizes.checksums[0] = (702381801U);
  helikopter_4_M->Sizes.checksums[1] = (673923214U);
  helikopter_4_M->Sizes.checksums[2] = (3704793831U);
  helikopter_4_M->Sizes.checksums[3] = (2608796625U);

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
