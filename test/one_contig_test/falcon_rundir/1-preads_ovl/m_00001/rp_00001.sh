export PATH=/nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/Falcon4.0/FALCON-integrate/fc_env/bin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/galaxy/home/szr165/bin:/galaxy/home/szr165/bin:${PATH}
export PYTHONPATH=:${PYTHONPATH}
export LD_LIBRARY_PATH=:/galaxy/home/szr165/R/usr/local/lib/R:/galaxy/home/szr165/R/usr/local/lib/R:${LD_LIBRARY_PATH}
set -vex
trap 'touch /nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/synth_human/try1/1-preads_ovl/m_00001/m_00001_done.exit' EXIT
cd /nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/synth_human/try1/1-preads_ovl/m_00001
hostname
date
time bash /nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/synth_human/try1/1-preads_ovl/m_00001/m_00001.sh
touch /nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/synth_human/try1/1-preads_ovl/m_00001/m_00001_done
