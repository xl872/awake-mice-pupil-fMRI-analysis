#!/bin/tcsh

# This script was created by the afni_proc.py quality control (APQC)
# generator.  
#
# It's purpose is to facilitate investigating the properties of the
# processing's residual dataset (errts*HEAD) file, by using the AFNI
# GUI's InstaCorr functionality.  

# As described in the popup help, users should just need to hold down
# the Ctrl+Shift keys and then left-click and move the mouse around
# (dragging or re-clicking).  Watch the correlation patterns to that
# seed location change, and this often provides an excellent way to
# understand the data.

# ver = 1.2
# -------------------------------------------------------------------------

set dset_ulay = "anat_final.allMice.average.stim+orig.HEAD"
set ic_dset   = "errts.allMice.average.stim+orig.HEAD"

set all_load  = ( "${dset_ulay}" "${ic_dset}"        \
                   "vr_base+orig.HEAD"               \
                   "mask_import_100mask+orig.HEAD"   \
                   all_runs.*.HEAD *.HEAD *.nii* )


set coord       = `3dinfo -dc3 "${ic_dset}"`

set voxvol      = `3dinfo -voxvol "${ic_dset}"`
set ic_seedrad  = `echo "${voxvol}"                                      \
                        | awk '{printf "%0.2f",(2*($1)^0.3334);}'`
echo "++ seedcorr radius: ${ic_seedrad}"

# ===========================================================================
# parameters set by default

setenv AFNI_THRESH_INIT_EXPON  0
setenv AFNI_NOSPLASH           YES
setenv AFNI_SPLASH_MELT        NO
setenv AFNI_STARTUP_WARNINGS   NO
setenv AFNI_NIFTI_TYPE_WARN    NO
setenv AFNI_NO_OBLIQUE_WARNING YES

# InstaCorr parameters

set ic_ignore   = 0
set ic_blur     = 0
set ic_automask = no
set ic_despike  = no
set ic_bandpass = 0,99999
set ic_polort   = -1
set ic_method   = P

# GUI visualization parameters

set pbar_sign   = "-"
set ncolors     = 99
set topval      = 0.6
set cbar        = "Reds_and_Blues_Inv"
set olay_alpha  = "Quadratic"
set olay_boxed  = "Yes"
set thresh      = 0.3
set frange      = ${topval}
set crossh      = MULTI
set xh_gap      = -1
set opacity     = 7
set OW          = "OPEN_WINDOW"

# port communication
set portnum = `afni -available_npb_quiet`

# ===========================================================================

afni -q  -no_detach                                                     \
    -npb ${portnum}                                                     \
     -com "SWITCH_UNDERLAY    ${dset_ulay}"                             \
     -com "INSTACORR INIT                                               \
                     DSET=${ic_dset}                                    \
                   IGNORE=${ic_ignore}                                  \
                     BLUR=${ic_blur}                                    \
                 AUTOMASK=${ic_automask}                                \
                  DESPIKE=${ic_despike}                                 \
                 BANDPASS=${ic_bandpass}                                \
                   POLORT=${ic_polort}                                  \
                  SEEDRAD=${ic_seedrad}                                 \
                   METHOD=${ic_method}"                                 \
     -com "INSTACORR SET      ${coord} J"                               \
     -com "SET_THRESHNEW      ${thresh}"                                \
     -com "SET_PBAR_ALL       ${pbar_sign}${ncolors} ${topval} ${cbar}" \
     -com "SET_FUNC_RANGE     ${frange}"                                \
     -com "SET_XHAIRS         ${crossh}"                                \
     -com "SET_XHAIR_GAP      ${xh_gap}"                                \
     -com "SET_FUNC_ALPHA     ${olay_alpha}"                            \
     -com "SET_FUNC_BOXED     ${olay_boxed}"                            \
     -com "$OW sagittalimage  opacity=${opacity}"                       \
     ${all_load:q}  &

sleep 1

set l = `prompt_user -pause \
"      Run InstaCorr on the residuals (errts) dataset\n\n\
\n\
InstaCorr calc using : ${ic_dset}\n\
Initial ulay dataset : ${dset_ulay}\n\
\n\
Wait briefly for the initial correlation patterns to appear.\n\
\n\
To use InstaCorr:\n\
Hold down Ctrl+Shift, and Left-click anywhere in the dataset.\n\
You can hold down Left-click and drag the cursor around, too.\n\
\n\
You will see the (unmasked) wholebrain correlation patterns\n\
from each clicked 'seed' location, updating instantly.\n\
Transparent thresholding is used (via 'A' and 'B' buttons).\n\
\n\
To jump to particular coordinates:\n\
+ Right-click -> 'Jump to (xyz)' \n\
+ Enter 3 space-separated coords\n\
+ Then, Right-click -> 'InstaCorr set',\n\
  or use standard Ctrl+Shift and Left-click.\n\
\n\
When done, hit 'OK' to exit.\n"`


if ("$l" != "1") then
    echo "+* Warn: InstaCorr guidance message failed to open"
endif

@Quiet_Talkers -npb_val ${portnum}

cat << EOF
===========================================
++ Goodbye, and thank you for InstaCorring.

EOF
exit 0

