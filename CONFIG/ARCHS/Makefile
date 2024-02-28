include ../Make.inc
arch=$(ARCH)
mach=$(ARCH)
machp=$(ARCH)p$(NPROC)
dcc=default
cc=gcc
skc=$(skcD)
dkc=$(dkcD)
smc=$(smcD)
dmc=$(dmcD)
comps = skc=$(skc) dkc=$(dkc) smc=$(smc) dmc=$(dmc)

adefd = $(ADuse)
anewd = $(BLDdir)/ARCHS/$(ADnew)

$(ADuse) : $(ADtar)
	bunzip2 -c $(ADtar) | tar xfm -
tarfile : $(ADnew)
	- rm $(ADnew).tar.bz2 $(ADnew).tar
	tar cvf $(ADnew).tar $(ADnew)
	$(BZIP) $(ADnew).tar

do_archdef : $(ADtarg)
do_nothing:
do_basic : basic.grd
do_full : full.grd
$(ADnew)/kern/$(dcc) :
	mkdir $(ADnew)/kern/$(dcc)
$(ADnew)/gemm/$(dcc) :
	mkdir $(ADnew)/gemm/$(dcc)
$(ADnew)/lapack/$(dcc) :
	mkdir $(ADnew)/lapack/$(dcc)

res:
	mkdir res
res/goodgcc.txt : res
	echo "'$(GOODGCC)'" > res/goodgcc.txt
$(ADnew) :
	mkdir $(ADnew)
	mkdir $(ADnew)/kern $(ADnew)/gemm $(ADnew)/lapack
	$(MAKE) $(ADnew)/kern/$(skc) dcc=$(skc)
	$(MAKE) $(ADnew)/kern/$(dkc) dcc=$(dkc)
	$(MAKE) $(ADnew)/gemm/$(smc) dcc=$(smc)
	$(MAKE) $(ADnew)/gemm/$(dmc) dcc=$(dmc)
	$(MAKE) $(ADnew)/lapack/$(dmc) dcc=$(dmc)
	$(MAKE) $(ADnew)/lapack/$(smc) dcc=$(smc)

xnegcp : $(SRCdir)/CONFIG/ARCHS/negcpfile.c
	$(XCC) $(XCCFLAGS) -o $@ $(SRCdir)/CONFIG/ARCHS/negcpfile.c
xnegr1 : $(SRCdir)/CONFIG/ARCHS/negr1file.c
	$(XCC) $(XCCFLAGS) -o $@ $(SRCdir)/CONFIG/ARCHS/negr1file.c
xnegmm : $(SRCdir)/CONFIG/ARCHS/negmmfile.c
	$(XCC) $(XCCFLAGS) -o $@ $(SRCdir)/CONFIG/ARCHS/negmmfile.c
xnegmv : $(SRCdir)/CONFIG/ARCHS/negmvfile.c
	$(XCC) $(XCCFLAGS) -o $@ $(SRCdir)/CONFIG/ARCHS/negmvfile.c
xnegflt : $(SRCdir)/CONFIG/ARCHS/negflt.c
	$(XCC) $(XCCFLAGS) -O -o xnegflt $(SRCdir)/CONFIG/ARCHS/negflt.c

ssys_put:
	- cp $(adefd)/kern/sMULADD $(SYSdir)/res/.
	- cp $(adefd)/kern/L1CacheSize $(SYSdir)/res/.
	- cp $(adefd)/kern/atlas_cache.h $(SYSdir)/res/.
	- cp $(adefd)/kern/sSIMD $(SYSdir)/res/.
ssys :
	- cp $(SYSdir)/res/atlas_cache.h $(anewd)/kern/.
	- cp $(SYSdir)/res/sMULADD $(anewd)/kern/.
	./xnegflt $(anewd)/kern/sMULADD
	- cp $(SYSdir)/res/L1CacheSize $(anewd)/kern/.
	- cp $(SYSdir)/res/sSIMD $(anewd)/res/.
skern_put:
	- cp $(adefd)/kern/$(skc)/[s,c]TRSM_NB $(BLDdir)/tune/blas/level3/res/.
	- cp $(adefd)/kern/$(skc)/[s,c]MV[N,T]K.sum $(BLDdir)/tune/blas/gemv/res/.
	- cp  $(adefd)/kern/$(skc)/[s,c]R?K.sum $(BLDdir)/tune/blas/ger/res/.
	- cp  $(adefd)/kern/$(skc)/[s,c]*.CPS $(BLDdir)/tune/blas/ger/res/.
	- cp $(adefd)/kern/$(skc)/[s,c]*_SUMM $(L1Tdir)/res/.
skern: xnegflt xnegr1 xnegmv xnegcp
	rm -f $(anewd)/kern/$(skc)/[s,c]*
	- cp $(BLDdir)/tune/blas/level3/res/[s,c]TRSM_NB $(anewd)/kern/$(skc)/.
	- cp $(BLDdir)/tune/blas/gemv/res/[s,c]MV[N,T]K.sum $(anewd)/kern/$(skc)/.
	- cp $(BLDdir)/tune/blas/ger/res/[s,c]R?K.sum $(anewd)/kern/$(skc)/.
	- cp $(BLDdir)/tune/blas/gemm/res/[s,c]cpylst.CPS $(anewd)/kern/$(skc)/.
	- cp $(BLDdir)/tune/blas/gemm/res/[s,c]cpyPERF.CPS $(anewd)/kern/$(skc)/.
	- ./xnegcp  $(anewd)/kern/$(skc)/[s,c]cpyPERF.CPS
	./xnegmv  $(anewd)/kern/$(skc)/[s,c]MVNK.sum
	./xnegmv  $(anewd)/kern/$(skc)/[s,c]MVTK.sum
	./xnegr1 $(anewd)/kern/$(skc)/[s,c]R?K.sum
	- cp $(L1Tdir)/res/[s,c]*_SUMM $(anewd)/kern/$(skc)/.
sgemmBAS_put:
	- cp $(adefd)/gemm/$(smc)/[s,c]CPYWORK.sum $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(smc)/[s,c]WORKING.sum $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(smc)/[s,c]AMMFRCLST.sum $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(smc)/[s,c]g[k,m]vAMMUR.sum $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(smc)/[s,c]opL1.sum $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(smc)/[s,c]gAMSYRK.sum $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(smc)/[s,c]gSYRKUM.sum $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(smc)/[s,c]AMSYRKFNL.sum $(MMTdir)/res/.
#	- cp $(adefd)/gemm/$(smc)/[s,c]ip*PERF.sum $(MMTdir)/res/.
sgemm_put:
	- cp $(adefd)/gemm/$(smc)/[s,c]*.sum $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(smc)/snreg $(MMTdir)/res/.
sgemm: xnegflt xnegmm
	rm -f $(anewd)/gemm/$(smc)/[s,c]*
	- cp $(MMTdir)/res/[s,c]trsm*.sum $(anewd)/gemm/$(smc)/.
	- cp $(MMTdir)/res/[s,c]trmm*.sum $(anewd)/gemm/$(smc)/.
	- cp $(MMTdir)/res/[s,c][i,o]p*.sum $(anewd)/gemm/$(smc)/.
	- cp $(MMTdir)/res/[s,c]*AMM*.sum $(anewd)/gemm/$(smc)/.
	- cp $(MMTdir)/res/[s,c]gSYRKUM.sum $(anewd)/gemm/$(smc)/.
	- cp $(MMTdir)/res/[s,c]gAMSYRK.sum $(anewd)/gemm/$(smc)/.
	- cp $(MMTdir)/res/sAMSYRKFNL.sum $(anewd)/gemm/$(smc)/.
	- cp $(MMTdir)/res/[s,c]SYRKFNL.sum $(anewd)/gemm/$(smc)/.
	- ./xnegmm $(anewd)/gemm/$(smc)/[s,c]*AMM*.sum \
                 $(anewd)/gemm/$(smc)/[s,c][i,o]p*.sum \
                 $(anewd)/gemm/$(smc)/[s,c]*SYRK*.sum 
	- cp $(MMTdir)/res/[s,c]WORKING.sum $(anewd)/gemm/$(smc)/.
	- cp $(MMTdir)/res/[s,c]CPYWORK.sum $(anewd)/gemm/$(smc)/.
	- cp $(MMTdir)/res/[s,c]nreg $(anewd)/gemm/$(smc)/.
	- cp $(MMTdir)/res/[s,c]fpuMM $(anewd)/gemm/$(smc)/.
sla_get:
	rm -f $(anewd)/lapack/$(smc)/atlas_[s,c]GetNB_*.h
	for file in $(INCAdir)/atlas_[s,c]GetNB_*.h ; do \
           if [ -s "$${file}" ]; then \
	      cp $${file} $(anewd)/lapack/$(smc)/. ; \
           fi \
        done
	rm -f $(anewd)/lapack/$(smc)/atlas_[s,c]tGetNB_*.h
	for file in $(INCAdir)/atlas_[s,c]tGetNB_*.h ; do \
           if [ -s "$${file}" ]; then \
	      cp $${file} $(anewd)/lapack/$(smc)/. ; \
           fi \
        done
sla_put :
	- cp $(adefd)/lapack/$(smc)/atlas_[s,c]GetNB_*.h $(INCAdir)/.
	- cp $(adefd)/lapack/$(smc)/atlas_[s,c]tGetNB_*.h $(INCAdir)/.
dsys_put:
	- cp $(adefd)/kern/dMULADD $(SYSdir)/res/.
	- cp $(adefd)/kern/L1CacheSize $(SYSdir)/res/.
	- cp $(adefd)/kern/atlas_cache.h $(SYSdir)/res/.
	- cp $(adefd)/kern/dSIMD $(SYSdir)/res/.
dsys :
	- cp $(SYSdir)/res/atlas_cache.h $(anewd)/kern/.
	- cp $(SYSdir)/res/dMULADD $(anewd)/kern/.
	./xnegflt $(anewd)/kern/dMULADD
	- cp $(SYSdir)/res/L1CacheSize $(anewd)/kern/.
	- cp $(SYSdir)/res/dSIMD $(anewd)/res/.
dkern_put:
	- cp $(adefd)/kern/$(dkc)/[d,z]TRSM_NB $(BLDdir)/tune/blas/level3/res/.
	- cp $(adefd)/kern/$(dkc)/[d,z]MV[N,T]K.sum $(BLDdir)/tune/blas/gemv/res/.
	- cp  $(adefd)/kern/$(dkc)/[d,z]R?K.sum $(BLDdir)/tune/blas/ger/res/.
	- cp  $(adefd)/kern/$(dkc)/[d,z]*.CPS $(BLDdir)/tune/blas/ger/res/.
	- cp $(adefd)/kern/$(dkc)/[d,z]*_SUMM $(L1Tdir)/res/.
dkern: xnegflt xnegr1 xnegmv xnegcp
	rm -f $(anewd)/kern/$(dkc)/[d,z]*
	- cp $(BLDdir)/tune/blas/level3/res/[d,z]TRSM_NB $(anewd)/kern/$(dkc)/.
	- cp $(BLDdir)/tune/blas/gemv/res/[d,z]MV[N,T]K.sum $(anewd)/kern/$(dkc)/.
	- cp $(BLDdir)/tune/blas/ger/res/[d,z]R?K.sum $(anewd)/kern/$(dkc)/.
	- cp $(BLDdir)/tune/blas/gemm/res/[d,z]cpylst.CPS $(anewd)/kern/$(dkc)/.
	- cp $(BLDdir)/tune/blas/gemm/res/[d,z]cpyPERF.CPS $(anewd)/kern/$(dkc)/.
	- ./xnegcp  $(anewd)/kern/$(dkc)/[d,z]cpyPERF.CPS
	./xnegmv  $(anewd)/kern/$(dkc)/[d,z]MVNK.sum
	./xnegmv  $(anewd)/kern/$(dkc)/[d,z]MVTK.sum
	./xnegr1 $(anewd)/kern/$(dkc)/[d,z]R?K.sum
	- cp $(L1Tdir)/res/[d,z]*_SUMM $(anewd)/kern/$(dkc)/.
dgemmBAS_put:
	- cp $(adefd)/gemm/$(dmc)/[d,z]CPYWORK.sum $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(dmc)/[d,z]WORKING.sum $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(dmc)/[d,z]AMMFRCLST.sum $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(dmc)/[d,z]g[k,m]vAMMUR.sum $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(dmc)/[d,z]opL1.sum $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(dmc)/[d,z]gAMSYRK.sum $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(dmc)/[d,z]gSYRKUM.sum $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(dmc)/[d,z]AMSYRKFNL.sum $(MMTdir)/res/.
#	- cp $(adefd)/gemm/$(dmc)/[d,z]ip*PERF.sum $(MMTdir)/res/.
dgemm_put:
	- cp $(adefd)/gemm/$(dmc)/[d,z]*.sum $(MMTdir)/res/.
	- cp $(adefd)/gemm/$(dmc)/dnreg $(MMTdir)/res/.
dgemm: xnegflt xnegmm
	rm -f $(anewd)/gemm/$(dmc)/[d,z]*
	- cp $(MMTdir)/res/[d,z]trsm*.sum $(anewd)/gemm/$(dmc)/.
	- cp $(MMTdir)/res/[d,z]trmm*.sum $(anewd)/gemm/$(dmc)/.
	- cp $(MMTdir)/res/[d,z][i,o]p*.sum $(anewd)/gemm/$(dmc)/.
	- cp $(MMTdir)/res/[d,z]*AMM*.sum $(anewd)/gemm/$(dmc)/.
	- cp $(MMTdir)/res/[d,z]gSYRKUM.sum $(anewd)/gemm/$(dmc)/.
	- cp $(MMTdir)/res/[d,z]gAMSYRK.sum $(anewd)/gemm/$(dmc)/.
	- cp $(MMTdir)/res/dAMSYRKFNL.sum $(anewd)/gemm/$(dmc)/.
	- cp $(MMTdir)/res/[d,z]SYRKFNL.sum $(anewd)/gemm/$(dmc)/.
	- ./xnegmm $(anewd)/gemm/$(dmc)/[d,z]*AMM*.sum \
                 $(anewd)/gemm/$(dmc)/[d,z][i,o]p*.sum \
                 $(anewd)/gemm/$(dmc)/[d,z]*SYRK*.sum 
	- cp $(MMTdir)/res/[d,z]WORKING.sum $(anewd)/gemm/$(dmc)/.
	- cp $(MMTdir)/res/[d,z]CPYWORK.sum $(anewd)/gemm/$(dmc)/.
	- cp $(MMTdir)/res/[d,z]nreg $(anewd)/gemm/$(dmc)/.
	- cp $(MMTdir)/res/[d,z]fpuMM $(anewd)/gemm/$(dmc)/.
dla_get:
	rm -f $(anewd)/lapack/$(dmc)/atlas_[d,z]GetNB_*.h
	for file in $(INCAdir)/atlas_[d,z]GetNB_*.h ; do \
           if [ -s "$${file}" ]; then \
	      cp $${file} $(anewd)/lapack/$(dmc)/. ; \
           fi \
        done
	rm -f $(anewd)/lapack/$(dmc)/atlas_[d,z]tGetNB_*.h
	for file in $(INCAdir)/atlas_[d,z]tGetNB_*.h ; do \
           if [ -s "$${file}" ]; then \
	      cp $${file} $(anewd)/lapack/$(dmc)/. ; \
           fi \
        done
dla_put :
	- cp $(adefd)/lapack/$(dmc)/atlas_[d,z]GetNB_*.h $(INCAdir)/.
	- cp $(adefd)/lapack/$(dmc)/atlas_[d,z]tGetNB_*.h $(INCAdir)/.

la_get : sla_get dla_get
la_put : sla_put dla_put
sys_put : ssys_put dsys_put
kern_put : skern_put dkern_put
gemm_put : sgemm_put dgemm_put
	- cp $(adefd)/gemm/atlas_cacheedge.h $(INCAdir)/.
	- cp $(adefd)/gemm/atlas_cacheedge.h $(MMTdir)/res/.
	- cp $(adefd)/gemm/atlas_tcacheedge.h $(INCAdir)/.
	- cp $(adefd)/gemm/atlas_tcacheedge.h $(MMTdir)/res/.
gemmBAS_put : sgemmBAS_put dgemmBAS_put
sys : ssys dsys
kern : skern dkern
gemm : sgemm dgemm
	- cp $(INCAdir)/atlas_cacheedge.h $(anewd)/gemm/
	- cp $(INCAdir)/atlas_tcacheedge.h $(anewd)/gemm/
ArchNew : $(ADnew) xnegflt
	$(MAKE) sys kern gemm la_get
	- cp $(BLDdir)/bin/INSTALL_LOG/?PerfSumm.txt $(anewd)/.
	rm -f xnegflt
basic.grd : $(ADuse)
	$(MAKE) sys_put kern_put gemmBAS_put
	touch basic.grd
	
full.grd : $(ADuse)
	$(MAKE) sys_put kern_put gemm_put la_put
	touch full.grd
archput : $(ADtarg)
