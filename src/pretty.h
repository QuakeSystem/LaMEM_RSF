/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : LaMEM
 **   License      : MIT, see LICENSE file for details
 **   Contributors : Anton Popov, Boris Kaus, see AUTHORS file for complete list
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : kaus@uni-mainz.de, popov@uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
//......................   PRETTY TERMINAL / RUN.LOG   ......................
//---------------------------------------------------------------------------
#ifndef __pretty_h__
#define __pretty_h__

struct FB;

//---------------------------------------------------------------------------

// Install PetscVFPrintf tee before PetscInitialize (always mirrors stdout/stderr to run.log once open).
void PrettyOutputInstallHook(void);

// Open/overwrite run.log next to the ParamFile; store case name from basename.
PetscErrorCode PrettyOutputOpenLog(void);

// Read pretty_output from .dat (default 0). If enabled, mute terminal verbose prints.
PetscErrorCode PrettyOutputReadOptions(FB *fb);

// Enter time-stepping solve stage (draw initial table; keep terminal muted if pretty).
void PrettyOutputEnterSolve(void);

// Leave solve stage (move cursor past the table).
void PrettyOutputLeaveSolve(void);

// Close run.log / restore printing.
PetscErrorCode PrettyOutputDestroy(void);

PetscInt PrettyOutputIsEnabled(void);
PetscInt PrettyOutputInSolve(void);

// Status setters (rank-0 table is redrawn when PrettyOutputRedraw is called).
void PrettyOutputSetStep(PetscInt step, PetscScalar time, PetscScalar dt, const char *time_unit);
void PrettyOutputSetDt(PetscScalar dt, const char *time_unit);
void PrettyOutputSetSNES(PetscInt snes_its, PetscScalar rel_res, PetscScalar fnorm, const char *reason, PetscInt ksp_its);
void PrettyOutputSetMRes(PetscScalar mRes2);
void PrettyOutputSetRSFControlled(PetscInt on);
void PrettyOutputSetCFLRestart(PetscInt on);
void PrettyOutputClearStepFlags(void);

// Redraw the fixed ANSI block (no-op unless pretty + solve stage + rank 0).
void PrettyOutputRedraw(void);

//---------------------------------------------------------------------------
#endif
