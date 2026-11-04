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
#include "LaMEM.h"
#include "parsing.h"
#include "pretty.h"
#include "tools.h"
#include <sys/ioctl.h>
#include <unistd.h>

//---------------------------------------------------------------------------

static const PetscInt PRETTY_WIDTH  = 160;
static const PetscInt PRETTY_NLINES = 9;

static FILE      *pretty_log       = NULL;
static PetscInt   pretty_enabled   = 0;
static PetscInt   pretty_mute      = 0;
static PetscInt   pretty_solve     = 0;
static PetscInt   pretty_drawn     = 0;
static char       pretty_case[_str_len_] = "LaMEM";
static PetscErrorCode (*pretty_vfprintf_default)(FILE *, const char[], va_list) = NULL;

static PetscInt   pretty_step      = 0;
static PetscScalar pretty_time     = 0.0;
static PetscScalar pretty_dt       = 0.0;
static char       pretty_tunit[_str_len_] = "s";
static PetscInt   pretty_snes_its  = 0;
static PetscScalar pretty_rel_res  = 0.0;
static PetscScalar pretty_fnorm    = 0.0;
static char       pretty_reason[_str_len_] = "-";
static PetscInt   pretty_ksp_its   = 0;
static PetscScalar pretty_mres2    = 0.0;
static PetscInt   pretty_rsf       = 0;
static PetscInt   pretty_cfl       = 0;

//---------------------------------------------------------------------------
static PetscInt PrettyTermCols(void)
{
	struct winsize ws;

	if(ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws) == 0 && ws.ws_col > 0)
	{
		return (PetscInt)ws.ws_col;
	}
	return PRETTY_WIDTH;
}
//---------------------------------------------------------------------------
// Keep one column free so "full-width line + \\n" does not create a blank row,
// and always emit an explicit newline (Cursor/IDE terminals often do not auto-wrap
// the same way ioctl reports).
static PetscInt PrettyDispWidth(PetscInt term_cols)
{
	PetscInt usable;

	if(term_cols < 2) usable = 1;
	else              usable = term_cols - 1;

	return (usable < PRETTY_WIDTH) ? usable : PRETTY_WIDTH;
}
//---------------------------------------------------------------------------
static void PrettyCropLine(char *line, PetscInt width)
{
	if((PetscInt)strlen(line) > width) line[width] = '\0';
}
//---------------------------------------------------------------------------
// One physical table row: content, erase-to-EOL (clears leftovers), newline.
static void PrettyEmitLine(FILE *out, const char *line)
{
	fputs(line, out);
	fputs("\033[K\n", out);
}
//---------------------------------------------------------------------------
static PetscErrorCode LaMEMVFPrintf(FILE *fd, const char format[], va_list Argp)
{
	char      buff[8192];
	size_t    length;
	FILE     *term;
	PetscBool is_stdout, is_stderr;

	PetscFunctionBeginUser;

	is_stdout = (PetscBool)(fd == stdout || (PETSC_STDOUT && fd == PETSC_STDOUT));
	is_stderr = (PetscBool)(fd == stderr || (PETSC_STDERR && fd == PETSC_STDERR));

	if(!is_stdout && !is_stderr)
	{
		// non-console PETSc viewers: pass through unchanged
		if(pretty_vfprintf_default) PetscCall((*pretty_vfprintf_default)(fd, format, Argp));
		PetscFunctionReturn(0);
	}

	PetscCall(PetscVSNPrintf(buff, sizeof(buff), format, &length, Argp));

	if(pretty_log)
	{
		fputs(buff, pretty_log);
		fflush(pretty_log);
	}

	// mute only normal stdout when pretty mode is on; keep stderr visible
	if(is_stderr || !pretty_mute)
	{
		term = is_stderr ? (PETSC_STDERR ? PETSC_STDERR : stderr) : (PETSC_STDOUT ? PETSC_STDOUT : stdout);
		fputs(buff, term);
		fflush(term);
	}

	PetscFunctionReturn(0);
}
//---------------------------------------------------------------------------
void PrettyOutputInstallHook(void)
{
	pretty_vfprintf_default = PetscVFPrintf;
	PetscVFPrintf = LaMEMVFPrintf;
}
//---------------------------------------------------------------------------
static void PrettyExtractCaseName(const char *paramfile, char *out)
{
	const char *base, *slash, *dot;
	size_t      n;

	slash = strrchr(paramfile, '/');
#ifdef _WIN32
	{
		const char *bslash = strrchr(paramfile, '\\');
		if(bslash && (!slash || bslash > slash)) slash = bslash;
	}
#endif
	base = slash ? slash + 1 : paramfile;

	dot = strrchr(base, '.');
	if(dot && dot != base) n = (size_t)(dot - base);
	else                   n = strlen(base);

	if(n >= (size_t)_str_len_) n = (size_t)_str_len_ - 1;
	memcpy(out, base, n);
	out[n] = '\0';
}
//---------------------------------------------------------------------------
static PetscErrorCode PrettyBuildLogPath(const char *paramfile, char *path)
{
	const char *slash;
	size_t      dirlen;

	PetscFunctionBeginUser;

	slash = strrchr(paramfile, '/');
#ifdef _WIN32
	{
		const char *bslash = strrchr(paramfile, '\\');
		if(bslash && (!slash || bslash > slash)) slash = bslash;
	}
#endif

	if(slash)
	{
		dirlen = (size_t)(slash - paramfile + 1);
		if(dirlen + strlen("run.log") >= (size_t)_str_len_)
		{
			SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "ParamFile path too long to form run.log path");
		}
		memcpy(path, paramfile, dirlen);
		path[dirlen] = '\0';
		strcat(path, "run.log");
	}
	else
	{
		PetscCall(PetscStrncpy(path, "run.log", (size_t)_str_len_));
	}

	PetscFunctionReturn(PETSC_SUCCESS);
}
//---------------------------------------------------------------------------
PetscErrorCode PrettyOutputOpenLog(void)
{
	PetscBool found;
	char      paramfile[_str_len_];
	char      logpath[_str_len_];

	PetscFunctionBeginUser;

	if(!ISRankZero(PETSC_COMM_WORLD)) PetscFunctionReturn(PETSC_SUCCESS);

	PetscCall(PetscOptionsGetCheckString("-ParamFile", paramfile, &found));
	if(found != PETSC_TRUE)
	{
		PetscFunctionReturn(PETSC_SUCCESS);
	}

	PrettyExtractCaseName(paramfile, pretty_case);
	PetscCall(PrettyBuildLogPath(paramfile, logpath));

	pretty_log = fopen(logpath, "w");
	if(!pretty_log)
	{
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Cannot open run.log at %s\n", logpath);
	}

	PetscFunctionReturn(PETSC_SUCCESS);
}
//---------------------------------------------------------------------------
PetscErrorCode PrettyOutputReadOptions(FB *fb)
{
	PetscInt val = 0;

	PetscFunctionBeginUser;

	PetscCall(getIntParam(fb, _OPTIONAL_, "pretty_output", &val, 1, 1));
	pretty_enabled = val;

	if(pretty_enabled && ISRankZero(PETSC_COMM_WORLD))
	{
		pretty_mute = 1;
	}

	PetscFunctionReturn(PETSC_SUCCESS);
}
//---------------------------------------------------------------------------
void PrettyOutputEnterSolve(void)
{
	pretty_solve = 1;
	pretty_drawn = 0;
	if(pretty_enabled) PrettyOutputRedraw();
}
//---------------------------------------------------------------------------
void PrettyOutputLeaveSolve(void)
{
	if(pretty_enabled && pretty_drawn && ISRankZero(PETSC_COMM_WORLD))
	{
		// leave the final table on screen; just stop overwriting
		fputc('\n', stdout);
		fflush(stdout);
	}
	pretty_solve = 0;
}
//---------------------------------------------------------------------------
PetscErrorCode PrettyOutputDestroy(void)
{
	PetscFunctionBeginUser;

	PrettyOutputLeaveSolve();

	if(pretty_log)
	{
		fclose(pretty_log);
		pretty_log = NULL;
	}

	pretty_mute    = 0;
	pretty_enabled = 0;
	pretty_solve   = 0;
	pretty_drawn   = 0;

	PetscFunctionReturn(PETSC_SUCCESS);
}
//---------------------------------------------------------------------------
PetscInt PrettyOutputIsEnabled(void)
{
	return pretty_enabled;
}
//---------------------------------------------------------------------------
PetscInt PrettyOutputInSolve(void)
{
	return pretty_solve;
}
//---------------------------------------------------------------------------
void PrettyOutputSetStep(PetscInt step, PetscScalar time, PetscScalar dt, const char *time_unit)
{
	pretty_step = step;
	pretty_time = time;
	pretty_dt   = dt;
	if(time_unit) { strncpy(pretty_tunit, time_unit, (size_t)_str_len_ - 1); pretty_tunit[_str_len_ - 1] = '\0'; }
}
//---------------------------------------------------------------------------
void PrettyOutputSetDt(PetscScalar dt, const char *time_unit)
{
	pretty_dt = dt;
	if(time_unit) { strncpy(pretty_tunit, time_unit, (size_t)_str_len_ - 1); pretty_tunit[_str_len_ - 1] = '\0'; }
}
//---------------------------------------------------------------------------
void PrettyOutputSetSNES(PetscInt snes_its, PetscScalar rel_res, PetscScalar fnorm, const char *reason, PetscInt ksp_its)
{
	pretty_snes_its = snes_its;
	pretty_rel_res  = rel_res;
	pretty_fnorm    = fnorm;
	pretty_ksp_its  = ksp_its;
	if(reason && reason[0]) { strncpy(pretty_reason, reason, (size_t)_str_len_ - 1); pretty_reason[_str_len_ - 1] = '\0'; }
	else                    { strncpy(pretty_reason, "-", (size_t)_str_len_ - 1); pretty_reason[_str_len_ - 1] = '\0'; }
}
//---------------------------------------------------------------------------
void PrettyOutputSetMRes(PetscScalar mRes2)
{
	pretty_mres2 = mRes2;
}
//---------------------------------------------------------------------------
void PrettyOutputSetRSFControlled(PetscInt on)
{
	pretty_rsf = on ? 1 : 0;
}
//---------------------------------------------------------------------------
void PrettyOutputSetCFLRestart(PetscInt on)
{
	pretty_cfl = on ? 1 : 0;
}
//---------------------------------------------------------------------------
void PrettyOutputClearStepFlags(void)
{
	// Only CFL is cleared on step accept. RSF stays sticky for the active dt and
	// is updated explicitly in ADVUpdateTimeStepRSF (avoids on/off blinking).
	pretty_cfl = 0;
}
//---------------------------------------------------------------------------
void PrettyOutputRedraw(void)
{
	char      line[PRETTY_WIDTH + 4];
	char      sep[PRETTY_WIDTH + 4];
	char      dt_extra[128];
	FILE     *out;
	PetscInt  term_cols, width;

	if(!pretty_enabled || !pretty_solve) return;
	if(!ISRankZero(PETSC_COMM_WORLD)) return;

	out       = PETSC_STDOUT ? PETSC_STDOUT : stdout;
	term_cols = PrettyTermCols();
	width     = PrettyDispWidth(term_cols);

	memset(sep, '-', (size_t)width);
	sep[width] = '\0';

	if(pretty_drawn)
	{
		// move to first row of previous block, column 0
		fprintf(out, "\033[%dA\r", PRETTY_NLINES);
	}

	dt_extra[0] = '\0';
	if(pretty_rsf) strcat(dt_extra, " (RSF controlled)");
	if(pretty_cfl) strcat(dt_extra, " (CFLMAX controlled restart)");

	snprintf(line, sizeof(line), "LaMEM  %s", pretty_case);
	PrettyCropLine(line, width);
	PrettyEmitLine(out, line);

	PrettyEmitLine(out, sep);

	snprintf(line, sizeof(line), "Step    %" PetscInt_FMT, pretty_step);
	PrettyCropLine(line, width);
	PrettyEmitLine(out, line);

	snprintf(line, sizeof(line), "Time    %.8e %s", (double)pretty_time, pretty_tunit);
	PrettyCropLine(line, width);
	PrettyEmitLine(out, line);

	snprintf(line, sizeof(line), "dt      %.4e %s%s", (double)pretty_dt, pretty_tunit, dt_extra);
	PrettyCropLine(line, width);
	PrettyEmitLine(out, line);

	snprintf(line, sizeof(line),
	         "SNES    %" PetscInt_FMT " iters, ||F||:%.4e, norm:%.4e (%s) %" PetscInt_FMT " KSP its",
	         pretty_snes_its, (double)pretty_rel_res, (double)pretty_fnorm, pretty_reason, pretty_ksp_its);
	PrettyCropLine(line, width);
	PrettyEmitLine(out, line);

	snprintf(line, sizeof(line), "|mRes|_2 %.4e", (double)pretty_mres2);
	PrettyCropLine(line, width);
	PrettyEmitLine(out, line);

	PrettyEmitLine(out, sep);

	snprintf(line, sizeof(line), "[full history -> run.log]");
	PrettyCropLine(line, width);
	PrettyEmitLine(out, line);

	fflush(out);
	pretty_drawn = 1;
}
//---------------------------------------------------------------------------
