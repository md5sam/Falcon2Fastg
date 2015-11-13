export PATH=/nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/Falcon4.0/FALCON-integrate/fc_env/bin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games:/galaxy/home/szr165/bin:/galaxy/home/szr165/bin:${PATH}
export PYTHONPATH=:${PYTHONPATH}
export LD_LIBRARY_PATH=:/galaxy/home/szr165/R/usr/local/lib/R:/galaxy/home/szr165/R/usr/local/lib/R:${LD_LIBRARY_PATH}
set -vex
trap 'touch /nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/synth_human/try1/2-asm-falcon/falcon_asm_done.exit' EXIT
cd /nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/synth_human/try1/1-preads_ovl
DB2Falcon -U preads
cd /nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/synth_human/try1/2-asm-falcon
find /nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/synth_human/try1/1-preads_ovl/las_files -name "*.las" > las.fofn 
fc_ovlp_filter --db /nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/synth_human/try1/1-preads_ovl/preads.db --fofn las.fofn --max_diff 1000 --max_cov 100 --min_cov 0 --bestn 1 --n_core 64 --min_len 100 > preads.ovl
ln -sf /nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/synth_human/try1/1-preads_ovl/preads4falcon.fasta .
fc_ovlp_to_graph preads.ovl --min_len 100 > fc_ovlp_to_graph.log
fc_graph_to_contig
fc_dedup_a_tigs
touch /nfs/brubeck.bx.psu.edu/scratch1/samarth/bonsai/synth_human/try1/2-asm-falcon/falcon_asm_done
