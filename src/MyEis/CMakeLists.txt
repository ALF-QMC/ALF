# Eispack library
SET(EIS_src ${SRCEIS}/balanc.f ${SRCEIS}/balbak.f ${SRCEIS}/cbabk2.f ${SRCEIS}/cbal.f ${SRCEIS}/cdiv.f ${SRCEIS}/cg.f ${SRCEIS}/ch.f
${SRCEIS}/comqr.f ${SRCEIS}/comqr2.f ${SRCEIS}/corth.f ${SRCEIS}/csroot.f ${SRCEIS}/elmhes.f ${SRCEIS}/eltran.f ${SRCEIS}/epslon.f ${SRCEIS}/hqr.f
 ${SRCEIS}/hqr2.f  ${SRCEIS}/htribk.f  ${SRCEIS}/htridi.f  ${SRCEIS}/pythag.f  ${SRCEIS}/rg.f  ${SRCEIS}/rs.f  ${SRCEIS}/tql2.f
  ${SRCEIS}/tqlrat.f  ${SRCEIS}/tred1.f  ${SRCEIS}/tred2.f)


ADD_LIBRARY(${MYEIS} STATIC ${EIS_src})
