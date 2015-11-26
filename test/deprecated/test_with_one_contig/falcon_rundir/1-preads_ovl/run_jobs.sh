# Daligner jobs (1)
daligner -v -h60 -t32 -H100 -e0.96 -l500 -s1000 preads.1 preads.1
# Initial sort jobs (1)
LAsort -v preads.1.preads.1.C0 preads.1.preads.1.N0 preads.1.preads.1.C1 preads.1.preads.1.N1 preads.1.preads.1.C2 preads.1.preads.1.N2 preads.1.preads.1.C3 preads.1.preads.1.N3 && LAmerge -v preads.1 preads.1.preads.1.C0.S preads.1.preads.1.N0.S preads.1.preads.1.C1.S preads.1.preads.1.N1.S preads.1.preads.1.C2.S preads.1.preads.1.N2.S preads.1.preads.1.C3.S preads.1.preads.1.N3.S && rm preads.1.preads.1.C0.S.las preads.1.preads.1.N0.S.las preads.1.preads.1.C1.S.las preads.1.preads.1.N1.S.las preads.1.preads.1.C2.S.las preads.1.preads.1.N2.S.las preads.1.preads.1.C3.S.las preads.1.preads.1.N3.S.las
