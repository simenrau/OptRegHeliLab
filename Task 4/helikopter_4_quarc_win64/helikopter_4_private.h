/*
 * helikopter_4_private.h
 *
 * Code generation for model "helikopter_4".
 *
<<<<<<< HEAD
 * Model version              : 1.185
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Tue Mar 06 17:40:03 2018
=======
 * Model version              : 1.180
 * Simulink Coder version : 8.6 (R2014a) 27-Dec-2013
 * C source code generated on : Tue Feb 27 21:35:31 2018
>>>>>>> 29b455799cc63d52fe3a4b40e62298412dcf4d14
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */
#ifndef RTW_HEADER_helikopter_4_private_h_
#define RTW_HEADER_helikopter_4_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"

/* Used by FromWorkspace Block: '<Root>/From Workspace' */
#ifndef rtInterpolate
# define rtInterpolate(v1,v2,f1,f2)    (((v1)==(v2))?((double)(v1)): (((f1)*((double)(v1)))+((f2)*((double)(v2)))))
#endif

#ifndef rtRound
# define rtRound(v)                    ( ((v) >= 0) ? floor((v) + 0.5) : ceil((v) - 0.5) )
#endif

#ifndef __RTWTYPES_H__
#error This file requires rtwtypes.h to be included
#endif                                 /* __RTWTYPES_H__ */

/* A global buffer for storing error messages (defined in quanser_common library) */
EXTERN char _rt_error_message[512];

/* private model entry point functions */
extern void helikopter_4_derivatives(void);

#endif                                 /* RTW_HEADER_helikopter_4_private_h_ */
