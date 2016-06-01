# Eispack library
SET(LIN_src ${SRCLIN}/cgedi.f
${SRCLIN}/cgefa.f
${SRCLIN}/dgedi.f
${SRCLIN}/dgefa.f
${SRCLIN}/zgedi.f
${SRCLIN}/zgefa.f
${SRCLIN}/zqrdc.f
${SRCLIN}/zqrsl.f
)

ADD_LIBRARY(${MYLIN} STATIC ${LIN_src})
