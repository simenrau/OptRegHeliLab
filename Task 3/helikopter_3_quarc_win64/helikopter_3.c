/*
 * helikopter_3.c
 *
 * Code generation for model "helikopter_3".
 *
 * Model version              : 1.177
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Tue Feb 20 14:13:13 2018
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
  real_T rtb_Sum4;
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

    /* Gain: '<S5>/Pitch: Count to rad' */
    helikopter_3_B.PitchCounttorad = helikopter_3_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S9>/Gain' */
    helikopter_3_B.Gain = helikopter_3_P.Gain_Gain *
      helikopter_3_B.PitchCounttorad;

    /* Gain: '<S5>/Travel: Count to rad' */
    helikopter_3_B.TravelCounttorad = helikopter_3_P.TravelCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o1;

    /* Gain: '<S12>/Gain' */
    helikopter_3_B.Gain_p = helikopter_3_P.Gain_Gain_a *
      helikopter_3_B.TravelCounttorad;

    /* SignalConversion: '<Root>/TmpSignal ConversionAtTo WorkspaceInport1' */
    helikopter_3_B.TmpSignalConversionAtToWorkspac[0] = helikopter_3_B.Gain;
    helikopter_3_B.TmpSignalConversionAtToWorkspac[1] = helikopter_3_B.Gain_p;
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
    /* Sum: '<Root>/Sum5' incorporates:
     *  Constant: '<Root>/Constant'
     */
    helikopter_3_B.Sum5 = helikopter_3_P.Constant_Value + helikopter_3_B.Gain_p;
  }

  /* TransferFcn: '<S5>/Travel: Transfer Fcn' */
  rtb_Backgain = 0.0;
  rtb_Backgain += helikopter_3_P.TravelTransferFcn_C *
    helikopter_3_X.TravelTransferFcn_CSTATE;
  rtb_Backgain += helikopter_3_P.TravelTransferFcn_D *
    helikopter_3_B.TravelCounttorad;

  /* Gain: '<S13>/Gain' */
  helikopter_3_B.Gain_d = helikopter_3_P.Gain_Gain_l * rtb_Backgain;

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
  rtb_Sum3_j[0] = helikopter_3_P.Gain1_Gain * helikopter_3_B.Sum5 - rtb_Sum3_j[0];
  rtb_Sum3_j[1] = helikopter_3_P.Gain1_Gain * helikopter_3_B.Gain_d -
    rtb_Sum3_j[1];
  rtb_Sum3_j[2] = helikopter_3_P.Gain1_Gain * helikopter_3_B.Gain - rtb_Sum3_j[2];
  rtb_Sum3_j[3] = helikopter_3_P.Gain1_Gain * helikopter_3_B.Gain_b -
    rtb_Sum3_j[3];

  /* Gain: '<Root>/Gain' */
  rtb_Backgain = ((helikopter_3_P.K[0] * rtb_Sum3_j[0] + helikopter_3_P.K[1] *
                   rtb_Sum3_j[1]) + helikopter_3_P.K[2] * rtb_Sum3_j[2]) +
    helikopter_3_P.K[3] * rtb_Sum3_j[3];

  /* Sum: '<Root>/Sum4' */
  rtb_Sum4 = rtb_Derivative - rtb_Backgain;
  if (rtmIsMajorTimeStep(helikopter_3_M)) {
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
  helikopter_3_B.Sum1 = ((rtb_Sum4 - helikopter_3_P.Gain1_Gain_f *
    helikopter_3_B.Gain) * helikopter_3_P.K_pp - helikopter_3_P.Gain1_Gain_f *
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
    rtb_Sum4 = helikopter_3_DW.TimeStampA;
    lastU = &helikopter_3_DW.LastUAtTimeA;
    if (helikopter_3_DW.TimeStampA < helikopter_3_DW.TimeStampB) {
      if (helikopter_3_DW.TimeStampB < helikopter_3_M->Timing.t[0]) {
        rtb_Sum4 = helikopter_3_DW.TimeStampB;
        lastU = &helikopter_3_DW.LastUAtTimeB;
      }
    } else {
      if (helikopter_3_DW.TimeStampA >= helikopter_3_M->Timing.t[0]) {
        rtb_Sum4 = helikopter_3_DW.TimeStampB;
        lastU = &helikopter_3_DW.LastUAtTimeB;
      }
    }

    rtb_Derivative = (helikopter_3_B.PitchCounttorad - *lastU) /
      (helikopter_3_M->Timing.t[0] - rtb_Sum4);
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
  rtb_Sum4 = (helikopter_3_B.Sum1 + helikopter_3_B.Sum2) *
    helikopter_3_P.Frontgain_Gain;

  /* Saturate: '<S5>/Front motor: Saturation' */
  if (rtb_Sum4 > helikopter_3_P.FrontmotorSaturation_UpperSat) {
    helikopter_3_B.FrontmotorSaturation =
      helikopter_3_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Sum4 < helikopter_3_P.FrontmotorSaturation_LowerSat) {
    helikopter_3_B.FrontmotorSaturation =
      helikopter_3_P.FrontmotorSaturation_LowerSat;
  } else {
    helikopter_3_B.FrontmotorSaturation = rtb_Sum4;
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
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.52359877559807444,
      0.52359877559803447, 0.52359877559798307, 0.52359877559791523,
      0.523598775597823, 0.5235987755976923, 0.52359877559749723,
      0.52359877559718215, 0.5235987755966055, 0.52359877559526735,
      0.52359877558925361, 0.38860546189573852, 0.10951698887526966,
      -0.11003833946102191, -0.27691123755291563, -0.397906827316469,
      -0.47963675127901129, -0.523591909441963, -0.52359877534342569,
      -0.52359877536343258, -0.5235987224145312, -0.503434893916266,
      -0.4649704404330644, -0.4207734586135255, -0.37347862555446437,
      -0.32522986550994049, -0.27772344682220723, -0.23225530276673972,
      -0.18977018817324479, -0.1509107408828376, -0.11606494657763664,
      -0.085410890684983476, -0.058958016011036271, -0.036584388943327453,
      -0.018069712651425214, -0.0031240162861197835, 0.008587900985857504,
      0.017426078578606414, 0.023758780967378854, 0.027949332008985622,
      0.030346034333594649, 0.031274796827770429, 0.031034093821064385,
      0.029891891654304288, 0.02808419891071845, 0.025814923189764541,
      0.023256747740391895, 0.020552773751695553, 0.017818707168261897,
      0.015145401405518497, 0.012601598411287872, 0.010236739517233149,
      0.0080837440173715247, 0.006161677142781272, 0.0044782499569109881,
      0.0030321116735972089, 0.0018149100878760932, 0.00081310836210016419,
      9.5565276842174477E-6, -0.00061517602288165625, -0.0010816948208518138,
      -0.0014108848952404623, -0.0016232063484833861, -0.0017381600630382566,
      -0.0017739020793310278, -0.0017469853873323276, -0.0016722086841820874,
      -0.0015625529081204719, -0.0014291879259837974, -0.0012815335105175721,
      -0.00112736059780635, -0.00097292068600022825, -0.00082309306297666382,
      -0.00068154128644839368, -0.00055087195231493544, -0.00043279025403042572,
      -0.00032824814515925642, -0.00023758206462032362, -0.00016063817112712341,
      -9.6883866645734262E-5, -4.5505078677972637E-5, -5.4893309923764264E-6,
      2.4304922656237182E-5, 4.5091880270843121E-5, 5.8113799538616E-5,
      6.4606933573472053E-5, 6.5776522273753557E-5, 6.2779201128075135E-5,
      5.6710998483779311E-5, 4.8598837981454217E-5, 3.9393103171466176E-5,
      2.995832301200616E-5, 2.1058409934228402E-5, 1.3332266149043414E-5,
      7.2554391803702265E-6, 3.0850995699961835E-6, 7.9190498794152707E-7,
      -4.0593743971645905E-17, -1.5332144927901074E-17, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0 } ;

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
      3.1415926535897931, 3.1378421413625279, 3.1262155534580049,
      3.1033093000299821, 3.0666274151912156, 3.014453922394225,
      2.9456562771176764, 2.8595077632937147, 2.7555515879654062,
      2.6335051104906531, 2.4931956060326312, 2.33451857606511,
      2.1583782719260483, 1.9678000778704734, 1.7674110445466211,
      1.5625810512761347, 1.3586839679386828, 1.1605920903248921,
      0.97235419679771284, 0.7969501232598406, 0.63643293374959053,
      0.49215959796927417, 0.36485729820996043, 0.25463993267206003,
      0.16109623593017916, 0.083400591252193068, 0.020423624533390569,
      -0.029167512815941696, -0.066823004653064913, -0.0940391853047939,
      -0.11230169021979586, -0.1230416759651909, -0.12760358353335186,
      -0.12722285379205736, -0.12301203948381058, -0.115953847476759,
      -0.10689976340781568, -0.096573044536239716, -0.085575006872611589,
      -0.074393673542786518, -0.0634139886026101, -0.05292893097858576,
      -0.043150984592951407, -0.034223531491125982, -0.026231833998759262,
      -0.019213359191285476, -0.013167274288011312, -0.0080630053551753542,
      -0.0038478045493674576, -0.00045331387359380022, 0.0021988530036483583,
      0.0041924638884134006, 0.0056128908698021623, 0.0065441252851624913,
      0.0070665002189294647, 0.0072550329904089855, 0.0071783005100552874,
      0.00689776339525422, 0.0064674596661119157, 0.00593399511017771,
      0.0053367645211473743, 0.0047083455725129364, 0.0040750147523033907,
      0.0034573422992259205, 0.0028708302460222221, 0.0023265643498006213,
      0.0018318567755375491, 0.0013908618413420552, 0.0010051519083085056,
      0.00067424460561774221, 0.00039607604483589989, 0.00016741753288223565,
      -1.5764411858470631E-5, -0.0001580031473864216, -0.00026407916771299981,
      -0.00033881565345682881, -0.00038691781724583337, -0.00041285121080772674,
      -0.00042075405089002686, -0.00041437869183475563, -0.00039705758179327461,
      -0.00037168935026994567, -0.000340741054982522, -0.00030626303876612115,
      -0.00026991328973278826, -0.00023298864146415995, -0.00019646057949526013,
      -0.00016101382346239438, -0.000127086220927492, -9.490881028218581E-5,
      -6.4545177863867948E-5, -3.5929439305404988E-5, -8.9023060705418331E-6,
      1.67552588483429E-5, 4.12913533646134E-5, 6.49561545002824E-5,
      8.7979679297890838E-5, 0.00011055569399486012, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.015002048909061993,
      -0.046506351618091038, -0.091625013712091113, -0.14672753935506733,
      -0.20869397118796096, -0.27519058110619465, -0.34459405529584691,
      -0.41582470131323507, -0.48818590989901217, -0.5612380178320886,
      -0.63470811987008435, -0.70456121655624571, -0.7623127762223002,
      -0.80155613329540887, -0.81931997308194582, -0.815588333349807,
      -0.79236751045516318, -0.7529515741087166, -0.70161629415148918,
      -0.642068758041, -0.57709334312126581, -0.50920919903725492,
      -0.44086946215160155, -0.3741747869675236, -0.31078257871194437,
      -0.25190786687521, -0.19836454939732906, -0.1506219673484929,
      -0.10886472260691592, -0.073050019660007862, -0.042959942981580158,
      -0.0182476302726438, 0.0015229189651780211, 0.016843257232987095,
      0.028232768028206252, 0.036216336275773309, 0.0413068754863039,
      0.043992150654512459, 0.044725333319300287, 0.043918739760705677,
      0.041940230496097367, 0.039111785542537421, 0.035709812407301693,
      0.03196678996946687, 0.028073899229895144, 0.024184339613096656,
      0.020417075731343829, 0.01686080322323159, 0.01357796270309463,
      0.010608667508968635, 0.0079744435390601675, 0.00568170792555505,
      0.0037249376614413143, 0.0020894997350678936, 0.00075413108591808494,
      -0.00030692992141479321, -0.001122148459204268, -0.0017212149165692209,
      -0.0021338582237368222, -0.0023889223561213389, -0.0025136757945377524,
      -0.0025333232808381796, -0.0024706898123098824, -0.0023460482128147939,
      -0.0021770635848864035, -0.00197883029705229, -0.0017639797367819752,
      -0.0015428397321341986, -0.0013236292107630541, -0.0011126742431273691,
      -0.000914634047814657, -0.00073272777896282509, -0.00056895494211180381,
      -0.00042430408130631294, -0.00029894594297531611, -0.0001924086551560182,
      -0.00010373357424757343, -3.1611360329200335E-5, 2.5501436221084996E-5,
      6.9284440165923942E-5, 0.00010147292609331565, 0.00012379318114969474,
      0.00013791206486560341, 0.00014539899613333143, 0.00014769859307451336,
      0.0001461122478755992, 0.00014178702413146298, 0.00013571041013960946,
      0.00012870964258122482, 0.00012145452967327142, 0.00011446295423385186,
      0.00010810853293945262, 0.00010263025967553892, 9.8144378065082008E-5,
      9.4659204542676026E-5, 9.20940991904337E-5, 9.0304058787877108E-5,
      8.9110958714119329E-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.10602875205861007, 0.22266037932307309, 0.31888147181624243,
      0.38944360631121488, 0.43795507377648218, 0.46997264230352059,
      0.49051724877498, 0.50343100141409236, 0.5114213858593829,
      0.51630439857559984, 0.51925862126752032, 0.49369500885904044,
      0.40816596705876168, 0.27735705984392112, 0.12554803518857244,
      -0.026373804426941448, -0.16411590764817854, -0.27857678423584409,
      -0.36281815260289824, -0.42085924264294028, -0.45922141703380015,
      -0.47977920385573458, -0.48299915977627905, -0.47137249195678615,
      -0.4480319169949194, -0.41610397764324641, -0.37842371849810491,
      -0.33742633592109816, -0.29512425777261331, -0.25312464196251705,
      -0.21266516986452463, -0.17465718802102484, -0.13973069118134346,
      -0.1082782996841416, -0.080496712382112964, -0.056424811263810687,
      -0.035977986944966435, -0.018978499319992549, -0.0051818550552452406,
      0.0057006952153475709, 0.013983347843646928, 0.019990368682992494,
      0.024043846827354082, 0.026454253043480879, 0.027513464961974845,
      0.027489921858446682, 0.02662558231072332, 0.025134375864333879,
      0.023201863003346496, 0.020985844389337645, 0.018617688948044889,
      0.016204179667314453, 0.013829704890780056, 0.011558650651908645,
      0.0094378756039130186, 0.00749917395600539, 0.00576165327421715,
      0.0042339729232725182, 0.0029164052970051357, 0.0018026958727821773,
      0.000881709656493849, 0.0001388609293404316, -0.00044266943438280026,
      -0.00088091922171180774, -0.0011943188110467248, -0.0014010371685183717,
      -0.0015184817035754474, -0.0015629330943505804, -0.0015492962434657843,
      -0.0014909491426511208, -0.0013996724643236189, -0.0012856440340190628,
      -0.0011574838622158488, -0.001022337038652467, -0.00088598344454119391,
      -0.00075296486124410322, -0.00062672160478955178, -0.00050973226282778129,
      -0.00040365143331496257, -0.00030944155013681613, -0.00022749592500501104,
      -0.00015775103811508149, -9.9786877774685778E-5, -5.2914770767496467E-5,
      -1.6252672910833923E-5, 1.1211682003878456E-5, 3.056902952022866E-5,
      4.2947186894513511E-5, 4.9478751347061E-5, 5.1276367137421249E-5,
      4.9413784960899316E-5, 4.4910622807850569E-5, 3.8718343149644832E-5,
      3.1704497960415815E-5, 2.4631875387714802E-5, 1.8129184956384727E-5,
      1.265132190729763E-5, 8.4323756487127923E-6, 5.4485883366234339E-6, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.42411500823444026,
      0.46652650905785203, 0.38488436997267733, 0.28224853797988991,
      0.1940458698610692, 0.12807027410815364, 0.082178425885837528,
      0.051655010556449575, 0.03196153778116162, 0.01953205086486786,
      0.011816890767682098, -0.10225444963391941, -0.34211616720111504,
      -0.523235628859362, -0.60723609862139483, -0.60768735846205546,
      -0.55096841288494836, -0.45784350635066229, -0.33696547346821648,
      -0.23216436016016828, -0.15344869756343951, -0.082231147287737713,
      -0.012879823682177783, 0.046506671277971562, 0.093362299847467037,
      0.12771175740669174, 0.15072103658056601, 0.16398953030802702,
      0.16920831259393937, 0.167998463240385, 0.16183788839196958,
      0.15203192737399915, 0.13970598735872553, 0.12580956598880746,
      0.11112634920811454, 0.096287604473209124, 0.081787297275377008,
      0.067997950499895543, 0.055186577058989245, 0.043530201082371246,
      0.033130610513197434, 0.024028083357382261, 0.016213912577446362,
      0.00964162486450719, 0.00423684767397586, -9.41724141126504E-5,
      -0.0034573581908934536, -0.00596482578555776, -0.0077300514439495272,
      -0.0088640744560354121, -0.009472621765171026, -0.0096540371229217464,
      -0.0094978991061376, -0.0090842169554856359, -0.0084831001919825084,
      -0.007754806591630514, -0.0069500827271529609, -0.0061107214037785292,
      -0.0052702705050695284, -0.0044548376968918347, -0.0036839448651533129,
      -0.0029713949086136696, -0.0023261214548929276, -0.0017529991493160301,
      -0.0012535983573396673, -0.000826873429886588, -0.00046977814022830267,
      -0.00017780556310053171, 5.4547403539183411E-5, 0.00023338840325865497,
      0.00036510671331000739, 0.00045611372121822372, 0.000512640687212856,
      0.00054058729425352761, 0.00054541437644509258, 0.00053207433318836275,
      0.00050497302581820568, 0.000467957367847082, 0.00042432331805127515,
      0.0003768395327125857, 0.00032778250052722041, 0.00027897954755971822,
      0.00023185664136158291, 0.00018748842802875725, 0.00014664839142665016,
      0.00010985741965884952, 7.7429390065400817E-5, 4.9512629497139395E-5,
      2.6126257810189975E-5, 7.1904631614409794E-6, -7.4503287060877088E-6,
      -1.8012648612194993E-5, -2.4769118632822942E-5, -2.805538075691609E-5,
      -2.8290490290804042E-5, -2.6010761725320305E-5, -2.1911452196348384E-5,
      -1.6875785034339347E-5, -1.1935149248357434E-5, -8.0237392002931063E-6,
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
  helikopter_3_M->Sizes.checksums[0] = (4170505535U);
  helikopter_3_M->Sizes.checksums[1] = (706050540U);
  helikopter_3_M->Sizes.checksums[2] = (45639923U);
  helikopter_3_M->Sizes.checksums[3] = (3958951230U);

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
    helikopter_3_B.PitchCounttorad = 0.0;
    helikopter_3_B.Gain = 0.0;
    helikopter_3_B.TravelCounttorad = 0.0;
    helikopter_3_B.Gain_p = 0.0;
    helikopter_3_B.TmpSignalConversionAtToWorkspac[0] = 0.0;
    helikopter_3_B.TmpSignalConversionAtToWorkspac[1] = 0.0;
    helikopter_3_B.Sum5 = 0.0;
    helikopter_3_B.Gain_d = 0.0;
    helikopter_3_B.Gain_b = 0.0;
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
  helikopter_3_M->Sizes.numBlocks = (63);/* Number of blocks */
  helikopter_3_M->Sizes.numBlockIO = (18);/* Number of block outputs */
  helikopter_3_M->Sizes.numBlockPrms = (147);/* Sum of parameter "widths" */
  return helikopter_3_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
