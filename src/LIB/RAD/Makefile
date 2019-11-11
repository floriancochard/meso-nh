OBJDIR=./$(ARCH)
MODDIR=$(ARCH)/MODULE_ECMWF_RAD2
VPATH=ECMWF_RAD:$(OBJDIR):.
BIBUSER=libecmwf_rad3b.a
INC = -I./ECMWF_RAD

# Juan : 24/11/2005 add objet "swclr_dst.o" for ECMWF_RAD2c
OBJS =  col2box.o legtri.o lwb.o lwbv.o lwc.o lw.o lwtt.o \
        lwttm.o lwu.o lwvb.o lwvd.o lwv.o lwvn.o olwb.o \
        olwbv.o olwc.o olw.o olwtt.o olwttm.o olwu.o olwvb.o \
        olwvd.o olwv.o olwvn.o parcli.o parclimf.o parcma.o \
        pardim.o pardimo.o parerob.o parfpos.o pargc5.o \
        parkind1.o parkind2.o parptrs.o parrint.o parrrtm.o \
        parrtm1d.o qsat.o rad1Driv_MACLATMOSPH_60LEVELS_ICRCCM3.o \
        radaca.o radaer.o radlsw.o radozc.o rrtm_cmbgb10.o \
        rrtm_cmbgb11.o rrtm_cmbgb12.o rrtm_cmbgb13.o rrtm_cmbgb14.o \
        rrtm_cmbgb15.o rrtm_cmbgb16.o rrtm_cmbgb1.o rrtm_cmbgb2.o \
        rrtm_cmbgb3.o rrtm_cmbgb4.o rrtm_cmbgb5.o rrtm_cmbgb6.o \
        rrtm_cmbgb7.o rrtm_cmbgb8.o rrtm_cmbgb9.o rrtm_ecrt_140gp.o \
        rrtm_gasabs1a_140gp.o rrtm_init_140gp.o rrtm_kgb10.o \
        rrtm_kgb11.o rrtm_kgb12.o rrtm_kgb13.o rrtm_kgb14.o \
        rrtm_kgb15.o rrtm_kgb16.o rrtm_kgb1.o rrtm_kgb2.o \
        rrtm_kgb3_a.o rrtm_kgb3_b.o rrtm_kgb3_c.o rrtm_kgb3.o \
        rrtm_kgb4_a.o rrtm_kgb4_b.o rrtm_kgb4_c.o rrtm_kgb4.o \
        rrtm_kgb5_a.o rrtm_kgb5_b.o rrtm_kgb5_c.o rrtm_kgb5.o \
        rrtm_kgb6.o rrtm_kgb7.o rrtm_kgb8.o rrtm_kgb9.o \
        rrtm_rrtm_140gp.o rrtm_rtrn1a_140gp.o rrtm_setcoef_140gp.o \
        rrtm_taumol10.o rrtm_taumol11.o rrtm_taumol12.o rrtm_taumol13.o \
        rrtm_taumol14.o rrtm_taumol15.o rrtm_taumol16.o rrtm_taumol1.o \
        rrtm_taumol2.o rrtm_taumol3.o rrtm_taumol4.o rrtm_taumol5.o \
        rrtm_taumol6.o rrtm_taumol7.o rrtm_taumol8.o rrtm_taumol9.o \
        satur.o strhandler.o suaerh.o suaerl.o suaersn.o suaerv.o \
        sucld.o sucldp.o suclop.o suclopn.o sucst.o suecaebc.o \
        suecaec.o suecaeor.o suecaesd.o suecaess.o suecaesu.o \
        suecozc.o suecozo.o sulw.o sulwn.o suolwn.o suovlp.o \
        surdi.o surrtab.o surrtftr.o surrtpk.o surrtrf.o susw.o \
        suswn.o sw1s.o sw2s.o swclr_dst.o swclr.o swde.o sw.o \
        swni.o swr.o swtt1.o swtt.o swu.o swuvo3.o yoeaerc.o \
        yoeaerd.o yoeaer.o yoeast.o yoecld.o yoecldp.o yoeclop.o \
        yoedbug.o yoelw.o yoemeth.o yoencst.o yoeolw.o yoeovlp.o \
        yoeozoc.o yoephli.o yoephy.o yoerad.o yoerdi.o yoerdu.o \
        yoerip.o yoerrta10.o yoerrta11.o yoerrta12.o yoerrta13.o \
        yoerrta14.o yoerrta15.o yoerrta16.o yoerrta1.o yoerrta2.o \
        yoerrta3.o yoerrta4.o yoerrta5.o yoerrta6.o yoerrta7.o \
        yoerrta8.o yoerrta9.o yoerrtab.o yoerrtam.o yoerrtbg2.o \
        yoerrtftr.o yoerrto10.o yoerrto11.o yoerrto12.o yoerrto13.o \
        yoerrto14.o yoerrto15.o yoerrto16.o yoerrto1.o yoerrto2.o \
        yoerrto3.o yoerrto4.o yoerrto5.o yoerrto6.o yoerrto7.o \
        yoerrto8.o yoerrto9.o yoerrtrf.o yoerrtrwt.o yoerrtwn.o \
        yoesat.o yoesw.o yoethf.o yoetluc.o yoevdf.o yoevdfs.o \
        yoeveg.o yoewcou.o yomcst.o yomlun.o yomrip.o

DIR_CONF:=$(shell pwd|sed -e 's/lib\/.*/conf/')

include $(DIR_CONF)/config.$(ARCH)
include Rules.$(ARCH)

$(BIBUSER): $(OBJDIR)/.dummy $(MODDIR)/.dummy $(OBJS)
	cd $(OBJDIR) ; pwd ; $(AR) r $@ $(OBJS) 

$(OBJDIR)/.dummy :
	-mkdir $(OBJDIR)
	@touch $(OBJDIR)/.dummy

$(MODDIR)/.dummy :
	-mkdir $(MODDIR)
	@touch $(MODDIR)/.dummy

clean:
	(if [ -d $(OBJDIR) ] ; then cd $(OBJDIR); rm -f *.f90 *.o ; fi)

distclean:
	-rm -rf $(BIBUSER) $(MODDIR) $(OBJDIR)

%.o:%.f90
	$(CPP) $(INC) $(CPPFLAGS) $< >  $(OBJDIR)/tmp.f90 
	$(SPL) $(OBJDIR)/tmp.f90 >  $(OBJDIR)/cpp_$*.f90 
	$(F90) $(INC) -I $(MODDIR) -c $(F90FLAGS) $(OBJDIR)/cpp_$*.f90
	rm -f tmp.f90 ; mv  cpp_$*.o $(OBJDIR)/$(*F).o
	-@mv  *.mod $(MODDIR)/. 2> /dev/null || echo pas de module dans $*.f90

# Juan : 24/11/2005 add rule "swclr_dst.o" for ECMWF_RAD2c
swclr_dst.o: swclr_dst.f90 parkind1.o \
	yoerad.o yoerdi.o yoerdu.o \
	yoesw.o 
# nombre de passe = 1
lw.o: lw.f90 parkind1.o \
	yoelw.o yoerdu.o 

lwb.o: lwb.f90 parkind1.o \
	yoelw.o 

lwbv.o: lwbv.f90 parkind1.o \
	yoelw.o yoerdu.o 

lwc.o: lwc.f90 parkind1.o \
	yoerad.o yoerdi.o yoeovlp.o

lwtt.o: lwtt.f90 parkind1.o \
	yoelw.o 

lwttm.o: lwttm.f90 parkind1.o \
	yoelw.o 

lwu.o: lwu.f90 parkind1.o \
	yoelw.o yoerdi.o yoerdu.o \
	yoesw.o yomcst.o 

lwv.o: lwv.f90 parkind1.o \
	yoelw.o 

lwvb.o: lwvb.f90 parkind1.o \
	yoelw.o 

lwvd.o: lwvd.f90 parkind1.o \
	yoelw.o 

lwvn.o: lwvn.f90 parkind1.o \
	yoelw.o 

parcli.o: parcli.f90 parkind1.o 

parclimf.o: parclimf.f90 parcli.o \
	parkind1.o 

parcma.o: parcma.f90 parkind1.o 

pardim.o: pardim.f90 parkind1.o 

pardimo.o: pardimo.f90 parkind1.o 

parerob.o: parerob.f90 parkind1.o 

parfpos.o: parfpos.f90 pardim.o \
	parkind1.o 

pargc5.o: pargc5.f90 parkind1.o 

parkind1.o: parkind1.f90 

parkind2.o: parkind2.f90 

parptrs.o: parptrs.f90 parkind1.o 

parrint.o: parrint.f90 parkind1.o 

parrrtm.o: parrrtm.f90 parkind1.o 

parrtm1d.o: parrtm1d.f90 parkind1.o 

radaer.o: radaer.f90 parkind1.o \
	yoeaerd.o 

radozc.o: radozc.f90 parkind1.o \
	yoeozoc.o yomcst.o 

rrtm_cmbgb1.o: rrtm_cmbgb1.f90 parkind1.o \
	parrrtm.o yoerrta1.o yoerrtftr.o \
	yoerrto1.o yoerrtrwt.o 

rrtm_cmbgb10.o: rrtm_cmbgb10.f90 parkind1.o \
	parrrtm.o yoerrta10.o yoerrtftr.o \
	yoerrto10.o yoerrtrwt.o 

rrtm_cmbgb11.o: rrtm_cmbgb11.f90 parkind1.o \
	parrrtm.o yoerrta11.o yoerrtftr.o \
	yoerrto11.o yoerrtrwt.o 

rrtm_cmbgb12.o: rrtm_cmbgb12.f90 parkind1.o \
	parrrtm.o yoerrta12.o yoerrtftr.o \
	yoerrto12.o yoerrtrwt.o 

rrtm_cmbgb13.o: rrtm_cmbgb13.f90 parkind1.o \
	parrrtm.o yoerrta13.o yoerrtftr.o \
	yoerrto13.o yoerrtrwt.o 

rrtm_cmbgb14.o: rrtm_cmbgb14.f90 parkind1.o \
	parrrtm.o yoerrta14.o yoerrtftr.o \
	yoerrto14.o yoerrtrwt.o 

rrtm_cmbgb15.o: rrtm_cmbgb15.f90 parkind1.o \
	parrrtm.o yoerrta15.o yoerrtftr.o \
	yoerrto15.o yoerrtrwt.o 

rrtm_cmbgb16.o: rrtm_cmbgb16.f90 parkind1.o \
	parrrtm.o yoerrta16.o yoerrtftr.o \
	yoerrto16.o yoerrtrwt.o 

rrtm_cmbgb2.o: rrtm_cmbgb2.f90 parkind1.o \
	parrrtm.o yoerrta2.o yoerrtftr.o \
	yoerrto2.o yoerrtrwt.o 

rrtm_cmbgb3.o: rrtm_cmbgb3.f90 parkind1.o \
	parrrtm.o yoerrta3.o yoerrtftr.o \
	yoerrto3.o yoerrtrwt.o 

rrtm_cmbgb4.o: rrtm_cmbgb4.f90 parkind1.o \
	parrrtm.o yoerrta4.o yoerrtftr.o \
	yoerrto4.o yoerrtrwt.o 

rrtm_cmbgb5.o: rrtm_cmbgb5.f90 parkind1.o \
	parrrtm.o yoerrta5.o yoerrtftr.o \
	yoerrto5.o yoerrtrwt.o 

rrtm_cmbgb6.o: rrtm_cmbgb6.f90 parkind1.o \
	parrrtm.o yoerrta6.o yoerrtftr.o \
	yoerrto6.o yoerrtrwt.o 

rrtm_cmbgb7.o: rrtm_cmbgb7.f90 parkind1.o \
	parrrtm.o yoerrta7.o yoerrtftr.o \
	yoerrto7.o yoerrtrwt.o 

rrtm_cmbgb8.o: rrtm_cmbgb8.f90 parkind1.o \
	parrrtm.o yoerrta8.o yoerrtftr.o \
	yoerrto8.o yoerrtrwt.o 

rrtm_cmbgb9.o: rrtm_cmbgb9.f90 parkind1.o \
	parrrtm.o yoerrta9.o yoerrtftr.o \
	yoerrto9.o yoerrtrwt.o 

rrtm_ecrt_140gp.o: rrtm_ecrt_140gp.f90 parkind1.o \
	parrrtm.o yoerad.o yoerdi.o \
	yoesw.o 

rrtm_gasabs1a_140gp.o: rrtm_gasabs1a_140gp.f90 parkind1.o \
	parrrtm.o yoerrtab.o 

rrtm_init_140gp.o: rrtm_init_140gp.f90 parkind1.o \
	parrrtm.o yoerrtbg2.o yoerrtftr.o \
	yoerrtrwt.o yoerrtwn.o 

rrtm_kgb1.o: rrtm_kgb1.f90 parkind1.o \
	yoerrto1.o 

rrtm_kgb10.o: rrtm_kgb10.f90 parkind1.o \
	yoerrto10.o 

rrtm_kgb11.o: rrtm_kgb11.f90 parkind1.o \
	yoerrto11.o 

rrtm_kgb12.o: rrtm_kgb12.f90 parkind1.o \
	yoerrta12.o yoerrto12.o 

rrtm_kgb13.o: rrtm_kgb13.f90 parkind1.o \
	yoerrta13.o yoerrto13.o 

rrtm_kgb14.o: rrtm_kgb14.f90 parkind1.o \
	yoerrto14.o 

rrtm_kgb15.o: rrtm_kgb15.f90 parkind1.o \
	yoerrta15.o yoerrto15.o 

rrtm_kgb16.o: rrtm_kgb16.f90 parkind1.o \
	yoerrta16.o yoerrto16.o 

rrtm_kgb2.o: rrtm_kgb2.f90 parkind1.o \
	yoerrta2.o yoerrto2.o 

rrtm_kgb6.o: rrtm_kgb6.f90 parkind1.o \
	yoerrto6.o 

rrtm_kgb7.o: rrtm_kgb7.f90 parkind1.o \
	yoerrta7.o yoerrto7.o 

rrtm_kgb8.o: rrtm_kgb8.f90 parkind1.o \
	yoerrta8.o yoerrto8.o 

rrtm_kgb9.o: rrtm_kgb9.f90 parkind1.o \
	yoerrta9.o yoerrto9.o 

rrtm_rrtm_140gp.o: rrtm_rrtm_140gp.f90 parkind1.o \
	parrrtm.o yoerrtwn.o 

rrtm_rtrn1a_140gp.o: rrtm_rtrn1a_140gp.f90 parkind1.o \
	parrrtm.o yoerrtab.o yoerrtftr.o \
	yoerrtwn.o 

rrtm_setcoef_140gp.o: rrtm_setcoef_140gp.f90 parkind1.o \
	parrrtm.o yoerrtrf.o 

rrtm_taumol1.o: rrtm_taumol1.f90 parkind1.o \
	parrrtm.o yoerrta1.o yoerrtwn.o 

rrtm_taumol10.o: rrtm_taumol10.f90 parkind1.o \
	parrrtm.o yoerrta10.o yoerrtwn.o 

rrtm_taumol11.o: rrtm_taumol11.f90 parkind1.o \
	parrrtm.o yoerrta11.o yoerrtwn.o 

rrtm_taumol12.o: rrtm_taumol12.f90 parkind1.o \
	parrrtm.o yoerrta12.o yoerrtwn.o 

rrtm_taumol13.o: rrtm_taumol13.f90 parkind1.o \
	parrrtm.o yoerrta13.o yoerrtwn.o 

rrtm_taumol14.o: rrtm_taumol14.f90 parkind1.o \
	parrrtm.o yoerrta14.o yoerrtwn.o 

rrtm_taumol15.o: rrtm_taumol15.f90 parkind1.o \
	parrrtm.o yoerrta15.o yoerrtwn.o 

rrtm_taumol16.o: rrtm_taumol16.f90 parkind1.o \
	parrrtm.o yoerrta16.o yoerrtwn.o 

rrtm_taumol2.o: rrtm_taumol2.f90 parkind1.o \
	parrrtm.o yoerrta2.o yoerrtbg2.o \
	yoerrtwn.o 

rrtm_taumol3.o: rrtm_taumol3.f90 parkind1.o \
	parrrtm.o yoerrta3.o yoerrtwn.o 

rrtm_taumol4.o: rrtm_taumol4.f90 parkind1.o \
	parrrtm.o yoerrta4.o yoerrtwn.o 

rrtm_taumol5.o: rrtm_taumol5.f90 parkind1.o \
	parrrtm.o yoerrta5.o yoerrtwn.o 

rrtm_taumol6.o: rrtm_taumol6.f90 parkind1.o \
	parrrtm.o yoerrta6.o yoerrtwn.o 

rrtm_taumol7.o: rrtm_taumol7.f90 parkind1.o \
	parrrtm.o yoerrta7.o yoerrtwn.o 

rrtm_taumol8.o: rrtm_taumol8.f90 parkind1.o \
	parrrtm.o yoerrta8.o yoerrtwn.o 

rrtm_taumol9.o: rrtm_taumol9.f90 parkind1.o \
	parrrtm.o yoerrta9.o yoerrtwn.o 

suaerh.o: suaerh.f90 parkind1.o \
	yoeaerd.o 

suaerl.o: suaerl.f90 parkind1.o \
	yoesw.o 

suaersn.o: suaersn.f90 parkind1.o \
	yoesw.o 

suaerv.o: suaerv.f90 parkind1.o 

sucld.o: sucld.f90 parkind1.o \
	yoecld.o 

sucldp.o: sucldp.f90 parkind1.o \
	yoecldp.o yomcst.o 

suclop.o: suclop.f90 parkind1.o \
	yoeclop.o 

suclopn.o: suclopn.f90 parkind1.o \
	yoesw.o 

sucst.o: sucst.f90 parkind1.o \
	yomcst.o yomrip.o 

suecaebc.o: suecaebc.f90 parkind1.o \
	yoeaerc.o 

suecaec.o: suecaec.f90 parkind1.o \
	yoeaerc.o yomcst.o 

suecaeor.o: suecaeor.f90 parkind1.o \
	yoeaerc.o 

suecaesd.o: suecaesd.f90 parkind1.o \
	yoeaerc.o 

suecaess.o: suecaess.f90 parkind1.o \
	yoeaerc.o 

suecaesu.o: suecaesu.f90 parkind1.o \
	yoeaerc.o 

suecozc.o: suecozc.f90 parkind1.o \
	yoeozoc.o yoerdi.o yomcst.o 

suecozo.o: suecozo.f90 parkind1.o \
	yoeozoc.o 

sulw.o: sulw.f90 parkind1.o \
	yoelw.o 

sulwn.o: sulwn.f90 parkind1.o \
	yoelw.o 

surdi.o: surdi.f90 parkind1.o \
	yoerdi.o 

surrtab.o: surrtab.f90 parkind1.o \
	yoerrtab.o 

surrtftr.o: surrtftr.f90 parkind1.o \
	yoerrtftr.o 

surrtpk.o: surrtpk.f90 parkind1.o \
	yoerrtwn.o 

surrtrf.o: surrtrf.f90 parkind1.o \
	yoerrtrf.o 

susw.o: susw.f90 parkind1.o \
	yoesw.o yomcst.o 

suswn.o: suswn.f90 parkind1.o \
	yoerad.o yoesw.o yomcst.o 

sw.o: sw.f90 parkind1.o \
	yoerad.o yoerdu.o 

swuvo3.o: swuvo3.f90 yoesw.o parkind1.o

sw1s.o: sw1s.f90 parkind1.o \
	yoerad.o yoesw.o 

sw2s.o: sw2s.f90 parkind1.o \
	yoerdu.o yoesw.o 

swclr.o: swclr.f90 parkind1.o \
	yoerad.o yoerdi.o yoerdu.o \
	yoesw.o 

swde.o: swde.f90 parkind1.o 

swni.o: swni.f90 parkind1.o \
	yoerad.o yoerdu.o yoesw.o 

swr.o: swr.f90 parkind1.o \
	yoecld.o yoerad.o yoeovlp.o 

swtt.o: swtt.f90 parkind1.o \
	yoesw.o 

swtt1.o: swtt1.f90 parkind1.o \
	yoesw.o 

swu.o: swu.f90 parkind1.o \
	yoecld.o yoerad.o yoerdu.o \
	yoesw.o yoeovlp.o

yoeaer.o: yoeaer.f90 parkind1.o 

yoeaerc.o: yoeaerc.f90 parkind1.o 

yoeaerd.o: yoeaerd.f90 parkind1.o 

yoeast.o: yoeast.f90 parkind1.o 

yoecld.o: yoecld.f90 parkind1.o 

yoecldp.o: yoecldp.f90 parkind1.o 

yoeclop.o: yoeclop.f90 parkind1.o 

yoedbug.o: yoedbug.f90 parkind1.o 

yoelw.o: yoelw.f90 parkind1.o 

yoemeth.o: yoemeth.f90 parkind1.o 

yoencst.o: yoencst.f90 parkind1.o 

yoeozoc.o: yoeozoc.f90 parkind1.o 

yoephli.o: yoephli.f90 parkind1.o 

yoephy.o: yoephy.f90 parkind1.o 

yoerad.o: yoerad.f90 parkind1.o 

yoerdi.o: yoerdi.f90 parkind1.o 

yoerdu.o: yoerdu.f90 parkind1.o 

yoerip.o: yoerip.f90 parkind1.o 

yoerrta1.o: yoerrta1.f90 parkind1.o 

yoerrta10.o: yoerrta10.f90 parkind1.o 

yoerrta11.o: yoerrta11.f90 parkind1.o 

yoerrta12.o: yoerrta12.f90 parkind1.o 

yoerrta13.o: yoerrta13.f90 parkind1.o 

yoerrta14.o: yoerrta14.f90 parkind1.o 

yoerrta15.o: yoerrta15.f90 parkind1.o 

yoerrta16.o: yoerrta16.f90 parkind1.o 

yoerrta2.o: yoerrta2.f90 parkind1.o 

yoerrta3.o: yoerrta3.f90 parkind1.o 

yoerrta4.o: yoerrta4.f90 parkind1.o 

yoerrta5.o: yoerrta5.f90 parkind1.o 

yoerrta6.o: yoerrta6.f90 parkind1.o 

yoerrta7.o: yoerrta7.f90 parkind1.o 

yoerrta8.o: yoerrta8.f90 parkind1.o 

yoerrta9.o: yoerrta9.f90 parkind1.o 

yoerrtab.o: yoerrtab.f90 parkind1.o 

yoerrtam.o: yoerrtam.f90 parkind1.o \
	parrrtm.o 

yoerrtbg2.o: yoerrtbg2.f90 parkind1.o \
	parrrtm.o 

yoerrtftr.o: yoerrtftr.f90 parkind1.o \
	parrrtm.o 

yoerrto1.o: yoerrto1.f90 parkind1.o 

yoerrto10.o: yoerrto10.f90 parkind1.o 

yoerrto11.o: yoerrto11.f90 parkind1.o 

yoerrto12.o: yoerrto12.f90 parkind1.o 

yoerrto13.o: yoerrto13.f90 parkind1.o 

yoerrto14.o: yoerrto14.f90 parkind1.o 

yoerrto15.o: yoerrto15.f90 parkind1.o 

yoerrto16.o: yoerrto16.f90 parkind1.o 

yoerrto2.o: yoerrto2.f90 parkind1.o 

yoerrto3.o: yoerrto3.f90 parkind1.o 

yoerrto4.o: yoerrto4.f90 parkind1.o 

yoerrto5.o: yoerrto5.f90 parkind1.o 

yoerrto6.o: yoerrto6.f90 parkind1.o 

yoerrto7.o: yoerrto7.f90 parkind1.o 

yoerrto8.o: yoerrto8.f90 parkind1.o 

yoerrto9.o: yoerrto9.f90 parkind1.o 

yoerrtrf.o: yoerrtrf.f90 parkind1.o 

yoerrtrwt.o: yoerrtrwt.f90 parkind1.o \
	parrrtm.o 

yoerrtwn.o: yoerrtwn.f90 parkind1.o 

yoesat.o: yoesat.f90 parkind1.o 

yoesw.o: yoesw.f90 parkind1.o 

yoethf.o: yoethf.f90 parkind1.o 

yoetluc.o: yoetluc.f90 parkind1.o 

yoevdf.o: yoevdf.f90 parkind1.o 

yoevdfs.o: yoevdfs.f90 parkind1.o 

yoeveg.o: yoeveg.f90 pardim.o \
	parkind1.o 

yoewcou.o: yoewcou.f90 parkind1.o 

yoeovlp.o: yoeovlp.f90 parkind1.o 

yomcst.o: yomcst.f90 parkind1.o 

yomlun.o: yomlun.f90 parkind1.o 

yomrip.o: yomrip.f90 parkind1.o 

rrtm_kgb3.o: rrtm_kgb3.f90 parkind1.o \
	yoerrta3.o yoerrto3.o 

rrtm_kgb3_a.o: rrtm_kgb3_a.f90 parkind1.o \
	yoerrta3.o yoerrto3.o 

rrtm_kgb3_b.o: rrtm_kgb3_b.f90 parkind1.o \
	yoerrta3.o yoerrto3.o 

rrtm_kgb3_c.o: rrtm_kgb3_c.f90 parkind1.o \
	yoerrta3.o yoerrto3.o 

rrtm_kgb4.o: rrtm_kgb4.f90 parkind1.o \
	yoerrta4.o yoerrto4.o 

rrtm_kgb4_a.o: rrtm_kgb4_a.f90 parkind1.o \
	yoerrta4.o yoerrto4.o 

rrtm_kgb4_b.o: rrtm_kgb4_b.f90 parkind1.o \
	yoerrta4.o yoerrto4.o 

rrtm_kgb4_c.o: rrtm_kgb4_c.f90 parkind1.o \
	yoerrta4.o yoerrto4.o 

rrtm_kgb5.o: rrtm_kgb5.f90 parkind1.o \
	yoerrta5.o yoerrto5.o 

rrtm_kgb5_a.o: rrtm_kgb5_a.f90 parkind1.o \
	yoerrta5.o yoerrto5.o 

rrtm_kgb5_b.o: rrtm_kgb5_b.f90 parkind1.o \
	yoerrta5.o yoerrto5.o 

rrtm_kgb5_c.o: rrtm_kgb5_c.f90 parkind1.o \
	yoerrta5.o yoerrto5.o 

