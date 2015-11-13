export PATH=/nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/Falcon4.0/FALCON-integrate/fc_env/bin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/galaxy/home/szr165/bin:/galaxy/home/szr165/bin:${PATH}
export PYTHONPATH=:${PYTHONPATH}
export LD_LIBRARY_PATH=:/galaxy/home/szr165/R/usr/local/lib/R:/galaxy/home/szr165/R/usr/local/lib/R:${LD_LIBRARY_PATH}
set -vex
trap 'touch /nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/synth_human/try1/1-preads_ovl/job_00000000/job_00000000_done.exit' EXIT
cd /nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/synth_human/try1/1-preads_ovl/job_00000000
hostname
date
time daligner_p -v -h60 -t32 -H100 -e0.96 -l500 -s1000 preads.1 preads.1
LAsort -v preads.1.preads.1.C0 preads.1.preads.1.N0 preads.1.preads.1.C1 preads.1.preads.1.N1 preads.1.preads.1.C2 preads.1.preads.1.N2 preads.1.preads.1.C3 preads.1.preads.1.N3 && LAmerge -v preads.1 preads.1.preads.1.C0.S preads.1.preads.1.N0.S preads.1.preads.1.C1.S preads.1.preads.1.N1.S preads.1.preads.1.C2.S preads.1.preads.1.N2.S preads.1.preads.1.C3.S preads.1.preads.1.N3.S && rm preads.1.preads.1.C0.S.las preads.1.preads.1.N0.S.las preads.1.preads.1.C1.S.las preads.1.preads.1.N1.S.las preads.1.preads.1.C2.S.las preads.1.preads.1.N2.S.las preads.1.preads.1.C3.S.las preads.1.preads.1.N3.S.las

rm -f preads.*.preads.*.*.las
 for f in `find $PWD -wholename "*.las"`; do mkdir -p ../m_00001; ln -sf $f ../m_00001; done 
touch /nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/synth_human/try1/1-preads_ovl/job_00000000/job_00000000_done
