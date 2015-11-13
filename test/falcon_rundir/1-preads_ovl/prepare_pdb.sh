export PATH=/nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/Falcon4.0/FALCON-integrate/fc_env/bin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/galaxy/home/szr165/bin:/galaxy/home/szr165/bin:${PATH}
export PYTHONPATH=:${PYTHONPATH}
export LD_LIBRARY_PATH=:/galaxy/home/szr165/R/usr/local/lib/R:/galaxy/home/szr165/R/usr/local/lib/R:${LD_LIBRARY_PATH}
set -vex
trap 'touch /nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/synth_human/try1/1-preads_ovl/pdb_build_done.exit' EXIT
cd /nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/synth_human/try1/1-preads_ovl
hostname
date
fasta2DB -v preads -f/nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/synth_human/try1/1-preads_ovl/preads.fofn
DBsplit -x100 -x500 -s50 preads
HPCdaligner -v -dal4 -t32 -h60 -e.96 -l500 -s1000 -H100 preads > /nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/synth_human/try1/1-preads_ovl/run_jobs.sh
touch /nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/synth_human/try1/1-preads_ovl/pdb_build_done
